function pathOptions(item, escape) {

  var markup = 
  "<div>" +
    "<div class = 'pull-left path-name' title = '" + escape(item.name) + "'>" +
        escape(item.name) +
    "</div>" +
    "<div class = 'pull-right path-fdr'>" +
        item.fdr +
    "</div>" +
    "<div class = 'clearfix'></div>" +
  "</div>";

  return markup;
}

