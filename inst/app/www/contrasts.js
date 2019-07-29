
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


function includeOptions(item, escape) {
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
  // styling if looking at cluster
  var res = "<div>" +
                     escape(item.label) +
                  "<span style='color: #A0A0A0;'>" +
                     " (" + item["pct.1"] + " :: " + item["pct.2"] + ")" +
                  "</span>" +
                "</div>";

  return res;

}
