$(document).on('shiny:sessioninitialized', function (e) {
  var mobile = /((iPhone)|(Android)|(BlackBerry))/.test(navigator.userAgent);
  Shiny.onInputChange('is_mobile', mobile);
});
