function pathOptions(item, escape) {

  var fdr = item.fdr ? item.fdr : '';

  var markup = 
  "<div>" +
    "<div class = 'pull-left path-name-option' title = '" + escape(item.name) + "'>" +
        escape(item.name) +
    "</div>" +
    "<div class = 'pull-right path-fdr'>" +
    fdr +
    "</div>" +
    "<div class = 'clearfix'></div>" +
  "</div>";

  return markup;
}

