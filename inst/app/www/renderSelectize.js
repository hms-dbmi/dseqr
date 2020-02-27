
// These functions are for styling the selectize input cluster selector

// styling if looking at cluster
function contrastOptions(item, escape) {

  var infoContrasts =  "<div style='color: #A0A0A0;text-align:right;'>" +
                        item.ncells + " :: " + item.pspace + item.pcells + "%" +
                       "</div>";

  var nsigf = typeof item.nsig == 'undefined' ? item.nbigf : item.nsigf;                        
  var infoSampleNsig = "<div style='color: #A0A0A0;text-align:right;'>" +
                        item.ntest + " :: " + item.nctrlf + " [<span style='color: dimgray;'>" + nsigf + "</span>]</div>";

  // for integrated dataset show number of cells in each test/ctrl sample
  var integratedTitle = item.ntest_each + " :: " + item.nctrl_each;
  var title = typeof item.ntest_each == 'undefined' ? '' : integratedTitle;

  var info = typeof item.ntest == 'undefined' ? infoContrasts : infoSampleNsig;

  var clustEl = "<div style='columns: 2;' title='" + title + "'>" +
                            "<div style='margin-right: -80px'>" +
                              "<div class='input-swatch' style='background-color:" + item.testColor + "'></div>" +
                              escape(item.name) +
                            "</div>" +
                            info +
                            "</div>";


  // styling if looking at contrast
  var conEl  = "<div title='"+ item.title + "'>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return (typeof item.pcells !== 'undefined' | typeof item.ntest !== 'undefined') ? clustEl : conEl;
}



//styling for current item
function contrastItem(item, escape) {
  var infoContrasts =  "<span style='color: #A0A0A0;'> (" + item.ncells + " :: " + item.pcells + "%" + ")</span>";

  var nsig = typeof item.nsig == 'undefined' ? item.nbig : item.nsig;
  var infoSampleNsig =  "<span style='color: #A0A0A0;'> (" + item.ntest + " :: " + item.nctrl + ") [<span style='color: dimgray;'>" + nsig + "</span>]</div>";
  var info = typeof item.ntest == 'undefined' ? infoContrasts : infoSampleNsig;

  // for integrated dataset show number of cells in each test/ctrl sample
  var integratedTitle = item.ntest_each + " :: " + item.nctrl_each;
  var title = typeof item.ntest_each == 'undefined' ? '' : integratedTitle;


  // styling if looking at cluster
  var clustEl = "<div title='" + title + "'>" +
                    "<div class='input-swatch' style='background-color:" + item.testColor + "'></div>" +
                     escape(item.name) +
                     info +
                "</div>";

   // styling if looking at contrast
  var conEl  = "<div title='"+ item.title + "'>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return (typeof item.pcells !== 'undefined' | typeof item.ntest !== 'undefined') ? clustEl : conEl;
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


function geneChoice(item, escape) {
  var res = "<div title = '" + item.description + "'>" + escape(item.label) + "</div>";
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

  var label = "<div class = 'path-fdr pull-left'>" + item.dirLabel + "</div>";
  var circle = item.ignore ? "<div></div>" : "<div class ='circle-swatch pull-right'></div>";  

  var markup = 
  "<div>" +
    "<div class = 'pull-left path-name-option' title = '" + escape(item.name) + "'>" +
        escape(item.label) +
    "</div>" +
    "<div class = 'pull-right path-option-right'>" +
      label +
      circle +
      "</div>" +
    "<div class = 'clearfix'></div>" +
  "</div>";

  return markup;
}

function pathItem(item, escape) {

  var fdr = item.fdr ? " (" + item.fdr + ")" : '';

  var markup = 
    "<div title = '" + escape(item.name) + "'>" +
        escape(item.label) +
        "<span class = 'path-fdr'>" + fdr + "</span>"
    "</div>";

  return markup;
}


function studyOption(item, escape) {

  var res = "<div style='columns: 2;'>" +
              "<div style='margin-right: -80px'>" +
                  escape(item.study) +
              "</div>" +
              "<div style='color: #A0A0A0;text-align:right;'>" + item.subset  + "</div>" +
            "</div>";


  return res;  
}


function studyItem(item, escape) {

  var res = "<div>" +
                  escape(item.study) +
                  "<span style='color: #A0A0A0;'>" + " (" + item.subset + ")" +"</span>"; +
            "</div>";

  return res;
}

function queryGenesItem(item, escape) {
  // color to indicate if in cmap alone or cmap and l1000
  var bgColor = item.cmap_only ? '#f0ad4e' : '#efefef';
  var res = "<div style='background-color:" + bgColor + "'>" + escape(item.gene) + "</div>";
  return res;
}

function queryGenesOption(item, escape) {

  var infoText = item.cmap_only ? 'ONLY CMAP02' : '';

  var res = "<div style='columns: 2;'>" +
              "<div style='margin-right: -80px'>" +
                  escape(item.gene) +
              "</div>" +
              "<div style='color: #A0A0A0;text-align:right;'>" + infoText + "</div>" +
            "</div>";

  return res;
}

function transferLabelOption(item, escape) {
  var predsMarkup = item.preds ? "<div class ='circle-swatch pull-right'></div>" : "<div></div>";

  var res = "<div style='columns: 2;'>" +
              "<div style='margin-right: -80px'>" +
                  escape(item.optionLabel) +
              "</div>" +
              predsMarkup +
            "</div>";

  return res;
}

function querySignatureItem(item, escape) {

  var info = item.type === "CMAP02/L1000 Perturbations" ? "" : " <span style='color: #A0A0A0;white-space: nowrap;'>" + "(" + item.type + ")" +"</span>";
 
  var res = "<div title='" + item.title + "'>" +
              escape(item.label) + info +
            "</div>";

  return res;
}

function querySignatureOptions(item, escape) {

  var res = "<div title='" + item.title + "'>" + escape(item.label) + "</div>";
  return res;
}


function pertOptions(item, escape) {

  var cor = item.cor ? item.cor : '';

  var markup = 
  "<div>" +
    "<div class = 'pull-left path-name-option'>" +
        escape(item.label) +
    "</div>" +
    "<div class = 'pull-right path-fdr'>" +
    cor+
    "</div>" +
    "<div class = 'clearfix'></div>" +
  "</div>";

  return markup;
}

function pertItem(item, escape) {

  var cor = item.cor ? item.cor : '';

  var markup = 
    "<div>" +
        escape(item.label) +
        "<span class = 'path-fdr'>" + " (" + cor + ")"  + "</span>" +
    "</div>";

  return markup;
}


function bulkContrastOptions(item, escape) {


  var clustEl = "<div>" +
                  "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
                  escape(item.name) +
                "</div>";

  return clustEl;
}



//styling for current item
function bulkContrastItem(item, escape) {

  var clustEl = "<div class='bulk-item'>" +
                  "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
                  escape(item.name) +
                  "<div class='contrast'>vs</div>" +
                "</div>";

  return clustEl;
}



function scDatasetOptions(item, escape) {

  var clustEl = "<div title='" + item.value + "'>" +
                  escape(item.optionLabel) +
                "</div>";

  return clustEl;
}



function scDatasetItem(item, escape) {

  var clustEl = "<div title='" + item.value + "'>" +
                  escape(item.itemLabel) +
                "</div>";

  return clustEl;
}
