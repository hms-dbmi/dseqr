$(document).on('DOMSubtreeModified', "#sc-form-dataset-up_raw_progress .progress-bar", function() {
    var value = $("#sc-form-dataset-up_raw_progress .progress-bar").html();
    Shiny.setInputValue('sc-form-dataset-up_raw_progress', value);
})

