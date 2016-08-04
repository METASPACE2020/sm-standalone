<!doctype HTML>
<html>
<head>
<title>IMS dataset browser</title>

<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css" integrity="sha384-fLW2N01lMqjakBkx3l/M9EahuwpSfeNvV63J5ezn3uZzapT0u7EYsXMjQV+0En5r" crossorigin="anonymous">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-slider/6.0.9/css/bootstrap-slider.min.css">

</head>
<body>
<div class="page-header container">

<form action="/show" method="get" class="text-center" role="form" class="form-horizontal">
<div class="col-sm-12">

  <div class="form-group col-sm-2" style="font-size:30px">
    {{!pretty_formula}}
  </div>

  <div class="form-group col-sm-10">
    <label class="col-sm-1">Tolerance (ppm):</label>
    <div class="col-sm-3">
      <input name="tolerance" type="text" id="tolerance" class="form-control" value="" data-slider-min="1" data-slider-max="10" data-slider-step="1" data-slider-value="{{tol}}" data-slider-orientation="horizontal" data-slider-selection="after" data-slider-tooltip="show">
    </div>

    <div class="col-sm-2">
      <input class="form-control" name="formula" type="text" value="{{formula}}" placeholder="Sum formula">
    </div>

    <div class="col-sm-4">
      <select class="form-control" name="dataset" required>
      % for dataset in datasets:
      <option value="{{dataset}}" {{"selected" if dataset == selected_dataset else ''}}>{{dataset}}</option>
      % end
      </select>
    </div>
  </div>

  <div class="form-group col-sm-2">
  </div>

  <div class="form-group col-sm-10">
    <label class="col-sm-1">Resolution:</label>
    <div class="col-sm-3">
      <input name="resolution" type="text" id="resolution" class="form-control" value="" data-slider-min="1000" data-slider-max="200000" data-slider-step="1000" data-slider-value="{{resolution}}" data-slider-orientation="horizontal" data-slider-selection="after" data-slider-tooltip="show">
    </div>

    <div class="col-sm-3">
      <input class="form-control btn-info" value="Show images" type="submit">
    </div>
    <div class="form-inline">
      <label class="col-sm-3 checkbox">
        <input class="form-control" type="checkbox" {{'checked="true"' if hs_removal else ''}} name="hs_removal">
        Remove hotspots
      </label>
    </div>
  </div>

  <div class="form-group col-sm-2">
  </div>

  <div class="form-group col-sm-10">
    <label class="col-sm-1">Number of peaks:</label>
    <div class="col-sm-3">
      <input name="npeaks" type="text" id="npeaks" class="form-control" value="" data-slider-min="4" data-slider-max="20" data-slider-step="1" data-slider-value="{{npeaks}}" data-slider-orientation="horizontal" data-slider-selection="after" data-slider-tooltip="show">
    </div>

    <div class="col-sm-6">
    </div>
  </div>
<!-- Pyisocalc cutoff: <input name="pyisocalc_cutoff" type="number" step="any" value="1e-5"><br/>-->
<!-- Points per FWHM: <input name="pts" type="number" step="1" value="10"></br> -->
    <input type="hidden" name="pyisocalc_cutoff" value="1e-3">
    <input type="hidden" name="pts" value="10">
    <input type="hidden" id="selected_adduct" name="adduct" value="{{selected_adduct}}">

  </div>

</div>
</form>

</div>

<div class="container">
  <ul class="nav nav-tabs">
  % for adduct, _ in isotope_patterns.iteritems():
    % if adduct == selected_adduct:
      <li class="nav active">
    % else:
      <li class="nav">
    % end

    <a data-toggle="tab" href="#adduct-{{adduct}}">+{{adduct}}</a></li>
  % end
  </ul>

  <div class="tab-content">
% for adduct, (mzs, intensities, stats) in isotope_patterns.iteritems():

  % if adduct == selected_adduct:
    <div id="adduct-{{adduct}}" class="tab-pane active">
  % else:
    <div id="adduct-{{adduct}}" class="tab-pane">
  % end

  <div class="text-center">
  <h4>
  % for score_name, score in stats.iteritems():
  {{score_name}}: {{"%.2f" % score}};
  % end
  </h4>
  </div>

  <table style='border-collapse:collapse'>
  <tr>
    <th>m/z</th>
    % for mz in mzs:
    <td style="border: 1px solid black" class="text-center">{{"%.4f" % mz}}</td>
    % end
  </tr>
  <tr>
    <th>intensity</th>
    % for intensity in intensities:
    <td style="border: 1px solid black" class="text-center">{{"%.1f" % intensity}}%</td>
    % end
  </tr>
  <tr>
    <th></th>
    % for mz in mzs:
    <td style='border: 1px solid black'>
    % if hs_removal:
      <img src="/show_image/{{selected_dataset}}/{{mz}}/{{tol}}?remove_hotspots=true" width="200px">
    % else:
      <img src="/show_image/{{selected_dataset}}/{{mz}}/{{tol}}" width="200px">
    % end
    </td>
    % end
  </tr>
  </table>
  <img src="/correlation_plot/{{selected_dataset}}/{{formula}}/{{adduct}}/{{",".join(map(str, mzs))}}/{{",".join(map(str, intensities))}}/{{tol}}"
       class="img-responsive" width="1100px"/>
  </div> <!-- adduct div -->

% end
  </div> <!-- tab content -->
</div> <!-- container -->

<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.0/jquery.min.js"></script>

<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-slider/6.0.9/bootstrap-slider.min.js"></script>
<script type="text/javascript">
  $("#tolerance").slider({formatter:function(ppm){return ''+ppm + " ppm";}});
  $("#resolution").slider({formatter:function(res){return ''+res;}});
  $("#npeaks").slider({formatter:function(n){return ''+n;}});

  $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
      $("#selected_adduct").val($(e.target).attr('href').split("-")[1]);
  });
</script>

<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>

</body>
</html>
