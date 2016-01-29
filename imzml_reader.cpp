#include <libxml2/libxml/xmlreader.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdint>

struct LibXmlString {
  xmlChar* str;

  LibXmlString(xmlChar* str) : str(str) {
  }

  bool operator==(const char* other) {
    return !xmlStrcmp(str, (xmlChar*)other);
  }

  bool operator==(const std::string& other) {
    return this->operator==(other.c_str());
  }

  bool operator!=(const char* other) {
    return !this->operator==(other);
  }

  bool operator!=(const std::string& other) {
    return !this->operator==(other);
  }

  operator std::string() const {
    return std::string((const char*)str);
  }

  ~LibXmlString() {
    xmlFree(str);
  }
};

LibXmlString getAttribute(xmlTextReaderPtr reader, const char* attribute) {
  return xmlTextReaderGetAttribute(reader, (xmlChar*)attribute);
}

class Metadata {
  static const std::map<std::string, std::string> supported_accessions_;

  std::map<std::string, std::string> dict_;
public:
  Metadata() {
  }

  void processNode(xmlTextReaderPtr reader) {
    auto accession = getAttribute(reader, "accession");
    if (accession == nullptr)
      return;

    auto it = supported_accessions_.find(accession);
    if (it == supported_accessions_.end())
      return;

    auto value = getAttribute(reader, "value");
    dict_[it->second] = value;
  }

  const std::map<std::string, std::string>& dict() {
    return dict_;
  }
};

const std::map<std::string, std::string> Metadata::supported_accessions_ =
{
  { "IMS:1000042", "max count of pixels x" },
  { "IMS:1000043", "max count of pixels y" },
  { "IMS:1000044", "max dimension x" },
  { "IMS:1000045", "max dimension y" },
  { "IMS:1000046", "pixel size x" },
  { "IMS:1000047", "pixel size y" },
  { "IMS:1000835", "matrix solution concentration" },
  {  "MS:1000843", "wavelength" },
  {  "MS:1000844", "focus diameter x" },
  {  "MS:1000845", "focus diameter y" },
  {  "MS:1000846", "pulse energy" },
  {  "MS:1000847", "pulse duration" },
  {  "MS:1000848", "attenuation" }
};

struct Spectrum {
  std::vector<double> mzs;
  std::vector<float> intensities;

  uint32_t x;
  uint32_t y;
  uint32_t z;
};

class ImzmlReader {
  Metadata metadata_;
  std::string filename_;
  xmlTextReaderPtr xml_reader_;
  int ret;

  std::string mz_group_name_, intensity_group_name_;

  void readMetadata() {
    std::string current_param_group;

    ret = xmlTextReaderRead(xml_reader_);
    while (ret == 1) {
      LibXmlString name = xmlTextReaderName(xml_reader_);
      if (name == "cvParam") {
        metadata_.processNode(xml_reader_);

        auto accession = getAttribute("accession");
        if (accession == "MS:1000514")
          mz_group_name_ = current_param_group;
        else if (accession == "MS:1000515")
          intensity_group_name_ = current_param_group;
      } else if (name == "spectrum" && isNodeStart()) {
        have_next_ = true;
        break;
      } else if (name == "referenceableParamGroup") {
        current_param_group = getAttribute("id");
      }

      ret = xmlTextReaderRead(xml_reader_);
    }
  }

  LibXmlString getAttribute(const char* attribute) {
    return ::getAttribute(xml_reader_, attribute);
  }

  template <typename T>
  T getIntValue() {
    auto value = getAttribute("value");
    return static_cast<T>(std::atoll((const char*)value.str));
  }

  template <typename T>
  void readIntValue(T& value) {
    value = getIntValue<T>();
  }

  struct ExternalArray {
    uint64_t file_offset;
    uint32_t length;
    uint32_t encoded_length;
  };

  ExternalArray mzs_;
  ExternalArray intensities_;

  bool have_next_;

  bool isNodeStart() const {
    return xmlTextReaderNodeType(xml_reader_) == XML_READER_TYPE_ELEMENT;
  }

  std::string ibd_filename_;

  template <typename T>
  void readExternalArray(const ExternalArray& array, std::vector<T>& buffer) {
    if (array.encoded_length / array.length == 4)
      readExternalArrayHelper<float, T>(array, buffer);
    else if (array.encoded_length / array.length == 8)
      readExternalArrayHelper<double, T>(array, buffer);
  }

  template <typename T, typename U>
  void readExternalArrayHelper(const ExternalArray& array, std::vector<U>& buffer) {
    std::ifstream in(ibd_filename_);
    in.seekg(array.file_offset);
    // FIXME: endianness?
    std::vector<T> vec(array.length);
    in.read(reinterpret_cast<char*>(&vec[0]), array.encoded_length);

    buffer.resize(array.length);
    std::copy(vec.begin(), vec.end(), buffer.begin());
  }

public:
  ImzmlReader(const std::string& filename) : filename_{filename}, have_next_{false} {
    xml_reader_ = xmlNewTextReaderFilename(filename.c_str());
    ibd_filename_ = filename_.substr(0, filename_.length() - 5) + "ibd";
    readMetadata();
  }

  bool readNextSpectrum(Spectrum& spectrum) {
    if (!have_next_)
      return false;

    have_next_ = false;

    spectrum.x = spectrum.y = spectrum.z = 0;
    spectrum.mzs.clear();
    spectrum.intensities.clear();

    ExternalArray* array;

    ret = xmlTextReaderRead(xml_reader_);
    while (ret == 1) {
      LibXmlString name = xmlTextReaderName(xml_reader_);

      if (name == "cvParam") {
        auto accession = getAttribute("accession");
        if      (accession == "IMS:1000050") readIntValue(spectrum.x);
        else if (accession == "IMS:1000051") readIntValue(spectrum.y);
        else if (accession == "IMS:1000052") readIntValue(spectrum.z);
        else if (accession == "IMS:1000102") readIntValue(array->file_offset);
        else if (accession == "IMS:1000103") readIntValue(array->length);
        else if (accession == "IMS:1000104") readIntValue(array->encoded_length);
      } else if (name == "spectrum" && isNodeStart()) {
        have_next_ = true;
        break;
      } else if (name == "referenceableParamGroupRef") {
        auto ref = getAttribute("ref");
        if (ref == mz_group_name_)
          array = &mzs_;
        else if (ref == intensity_group_name_)
          array = &intensities_;
      }
      ret = xmlTextReaderRead(xml_reader_);
    }

    readExternalArray<double>(mzs_, spectrum.mzs);
    readExternalArray<float>(intensities_, spectrum.intensities);

    return true;
  }

  ~ImzmlReader() {
    xmlFreeTextReader(xml_reader_);
  }

  const std::map<std::string, std::string>& dict() {
    return metadata_.dict();
  }
};
