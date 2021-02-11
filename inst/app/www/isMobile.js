$(document).on('shiny:sessioninitialized', function (e) {
  var mobile = /((iPhone)|(iPad)|(Android)|(BlackBerry))/.test(navigator.userAgent);
  Shiny.onInputChange('is_mobile', mobile);
});
