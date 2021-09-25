
// These functions are for styling the selectize input cluster selector
// Item: selected item
// Options: available options in dropdown


// Single Cell

// selected cluster or contrast
function contrastOptions(item, escape) {
    
    var infoContrasts =
    "<div class='selectize-info'>" +
    "%" +
    item.pcells +
    " :: " +
    item.ncellsf +
    "</div>";
    
    var nsigf = typeof item.nsig == 'undefined' ? item.nbigf : item.nsigf;
    
    var infoSampleNsig = 
    "<div class='selectize-info'>" +
    "<span>[</span>" + 
    "<span style='color: dimgray;'>" + nsigf + "</span>" +
    "<span>]</span> " + 
    item.nctrlf + 
    "<span> :: </span>" +
    item.ntest +
    "</div>";
    
    // for integrated dataset show number of cells in each test/ctrl sample as title element
    // otherwise show potentially truncated label
    var integratedTitle = item.ntest_each + " :: " + item.nctrl_each;
    var title = typeof item.ntest_each == 'undefined' ? item.label : integratedTitle;
    
    var info = typeof item.ntest == 'undefined' ? infoContrasts : infoSampleNsig;
    
    var swatchClass = item.testColor == '' ? '' : 'input-swatch';
    
    // disable when no top_table
    var disabled = item.disabled ? 'disabled-option': '';
    
    var clustEl = 
    "<div style='columns: 2;' title='" + title + "' class='" + disabled + "'>" +
    "<div style='margin-right: -80px'>" +
    "<div class='" + swatchClass +"' style='background-color:" + item.testColor + "'></div>" +
    escape(item.name) +
    "</div>" + info + "</div>";
    
    
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

function contrastItem(item, escape) {
    var infoContrasts =  "<span style='color: #A0A0A0;'> (" + item.ncells + " :: " + item.pcells + "%" + ")</span>";
    
    var nsig = typeof item.nsig == 'undefined' ? item.nbig : item.nsig;
    var infoSampleNsig =  "<span style='color: #A0A0A0;'> (" + item.ntest + " :: " + item.nctrl + ") [<span style='color: dimgray;'>" + nsig + "</span>]</div>";
    var info = typeof item.ntest == 'undefined' ? infoContrasts : infoSampleNsig;
    
    // for integrated dataset show number of cells in each test/ctrl sample
    var integratedTitle = item.ntest_each + " :: " + item.nctrl_each;
    var title = typeof item.ntest_each == 'undefined' ? '' : integratedTitle;
    
    console.log(item.testColor == '');
    
    var swatchClass = item.testColor == '' ? '' : 'input-swatch';
    
    // disable when no top_table
    var disabled = item.disabled ? 'disabled-option': '';
    
    // styling if looking at cluster
    var clustEl = "<div title='" + title + "' class='" + disabled + "'>" +
    "<div class='" + swatchClass +"' style='background-color:" + item.testColor + "'></div>" +
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


function transferLabelOption(item, escape) {
    var predsMarkup = item.preds ? "<div class ='circle-swatch pull-right'></div>" : "<div></div>";
    
    var res = "<div title='" + item.label + "' style='columns: 2;'>" +
    "<div style='margin-right: -80px'>" +
    escape(item.optionLabel) +
    "</div>" +
    predsMarkup +
    "</div>";
    
    return res;
}


function scDatasetOptions(item, escape) {
    
    var opt = item.optgroup;
    if (opt === item.label) {
        opt = ''
    } else {
        opt = typeof opt == 'undefined' ? '' : opt;
        opt = (opt == 'Previous Session' || opt == 'Integrated' || opt =='Individual') ? '' : opt + '_';
    }
    
    var full_name = opt + item.label;
    
    var clustEl = "<div title='" + full_name + "'>" +
    escape(item.label) +
    "</div>";
    
    return clustEl;
}


function scDatasetItem(item, escape) {
    
    
    var opt = item.optgroup;
    if (opt === item.label) {
        opt = ''
    } else {
        opt = typeof opt == 'undefined' ? '' : opt;
        opt = (opt == 'Previous Session' || opt == 'Integrated' || opt =='Individual') ? '' : opt + '_';
    }
    
    var full_name = opt + item.label;
    
    var clustEl = "<div title='" + full_name + "'>" +
    escape(full_name) +
    "</div>";
    
    return clustEl;
}


function scDatasetItemDF(item, escape) {
    
    
    var label = item.label;
    
    var clustEl = "<div title='" + label + "'>" +
    escape(label) +
    "</div>";
    
    return clustEl;
}



// Drugs tab:

// select cell line
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


function querySignatureItem(item, escape) {
    
    var info = item.group === "CMAP02/L1000 Perturbations" ? "" : " <span style='color: #A0A0A0;white-space: nowrap;'>" + "(" + item.group + ")" +"</span>";
    
    var res = "<div title='" + item.title + "'>" +
    escape(item.dataset_name) + info +
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


// Bulk tab:
function bulkContrastOptions(item, escape) {
    
    var swatch = item.color ? "<div class='input-swatch' style='background-color:" + item.color + "'></div>" : "";
    
    var clustEl = "<div>" +
    swatch +
    escape(item.name) +
    "</div>";
    
    return clustEl;
}


function bulkContrastItem(item, escape) {
    
    var swatch = "<div class='input-swatch' style='background-color:" + item.color + "'></div>";
    
    var clustEl = "<div class='bulk-item'>" +
    swatch +
    escape(item.name) +
    "<div class='contrast'>vs</div>" +
    "</div>";
    
    return clustEl;
}


function bulkDatasetOptions(item, escape) {
    
    var clustEl = "<div title='" + item.label + "'>" +
    escape(item.value) +
    "</div>";
    
    return clustEl;
}


function bulkDatasetItem(item, escape) {
    
    var clustEl = "<div>" +
    escape(item.value) +
    "</div>";
    
    return clustEl;
}
