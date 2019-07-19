// styling if looking at cluster
function cellOptions(item, escape) {

  var markup = "<div style='columns: 2;'>" +
    "<div>" +
        escape(item.cell_id) +
    "</div>" +
    "<div style='color: #A0A0A0;text-align:right;'>" +
        item.sample_type +
    "</div>" +
  "</div>";

  return markup;
}

