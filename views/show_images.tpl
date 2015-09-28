<h2>{{formula}}</h2>
% for adduct, (mzs, intensities, stats) in isotope_patterns.iteritems():
  <h3>+{{adduct}}</h3>
  % for score_name, score in stats.iteritems():
  {{score_name}}: {{score}};
  % end
  <table style='border-collapse:collapse'>
  <tr>
    <th>m/z</th>
    % for mz in mzs:
    <td style="border: 1px solid black">{{mz}}</td>
    % end
  </tr>
  <tr>
    <th>intensity</th>
    % for intensity in intensities:
    <td style="border: 1px solid black">{{intensity}}</td>
    % end
  </tr>
  <tr>
    <th></th>
    % for mz in mzs:
    <td style='border: 1px solid black'>
    % if hs_removal:
      <img src="/show_image/{{mz}}/{{tol}}?remove_hotspots=true" width="200px">
    % else:
      <img src="/show_image/{{mz}}/{{tol}}" width="200px">
    % end
    </td>
    % end
  </tr>
  </table>
% end
