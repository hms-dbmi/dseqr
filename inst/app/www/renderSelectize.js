
// These functions are for styling the selectize input cluster selector

// styling if looking at cluster
function contrastOptions(item, escape) {

  var clustEl = "<div style='columns: 2;'>" +
                  "<div style='margin-right: -80px'>" +
                    "<div class='input-swatch' style='background-color:" + item.testColor + "'></div>" +
                     escape(item.name) +
                  "</div>" +
                  "<div style='color: #A0A0A0;text-align:right;'>" +
                   item.ncells + " :: " + item.pspace + item.pcells + "%" +
                  "</div>" +
                "</div>";

  // styling if looking at contrast
  var conEl  = "<div>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return typeof item.pcells !== 'undefined' ? clustEl : conEl;
}



//styling for current item
function contrastItem(item, escape) {
  // styling if looking at cluster
  var clustEl = "<div>" +
                    "<div class='input-swatch' style='background-color:" + item.testColor + "'></div>" +
                     escape(item.name) +
                  "<span style='color: #A0A0A0;'>" +
                     " (" + item.ncells + " :: " + item.pcells + "%)" +
                  "</span>" +
                "</div>";

   // styling if looking at contrast
  var conEl  = "<div>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return typeof item.pcells !== 'undefined' ? clustEl : conEl;
}


function excludeOptions(item, escape) {
  var res = "<div>" +
                "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
                  escape(item.name) +
            "</div>";

  return res;
}

function integationOption(item, escape) {
  var res = "<div>" +
              "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
                escape(item.label) +
            "</div>";

  return res;
}

function geneOption(item, escape) {
  var res = "<div style='columns: 2;'>" +
                  "<div style='margin-right: -80px'>" +
                     escape(item.label) +
                  "</div>" +
                  "<div style='color: #A0A0A0;text-align:right;'>" +
                   item["pct.1"] + " :: " + item.pspace + item["pct.2"] +
                  "</div>" +
                "</div>";


  return res;
}


function geneItem(item, escape) {
  var res = "<div>" +
                     escape(item.label) +
                  "<span style='color: #A0A0A0;'>" +
                     " (" + item["pct.1"] + " :: " + item["pct.2"] + ")" +
                  "</span>" +
                "</div>";

  return res;

}

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




function studyOption(item, escape) {

  subsetText = item.subset == null ? "" : item.subset;

  var res = "<div style='columns: 2;'>" +
                  "<div style='margin-right: -80px'>" +
                     escape(item.study) +
                  "</div>" +
                  "<div style='color: #A0A0A0;text-align:right;'>" +
                    subsetText  +
                  "</div>" +
                "</div>";


  return res;  
}


function studyItem(item, escape) {
  console.log(item);

  var subsetMarkup = item.subset == null ? 
  "" : 
  "<span style='color: #A0A0A0;'>" + " (" + item.subset + ")" +"</span>";

  var res = "<div>" +
                  escape(item.study) +
                  subsetMarkup +
            "</div>";

  return res;
}