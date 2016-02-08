void* imzb_reader_new(const char*);
void imzb_reader_free(void*);
int imzb_reader_height(void*);
int imzb_reader_width(void*);
void imzb_reader_image(void*, double, double, float*);
