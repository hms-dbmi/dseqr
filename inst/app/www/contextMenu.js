var title = null;
var menuHeight = 0;

function setupContextMenu() {
  title = null;
  menuHeight = $("#cor-menu").height();

  $('.cor-point').off().on("contextmenu", function(event) {
    // Avoid the real one
    event.preventDefault();

    
    corpoint = $(event.target).closest('.cor-point');
    corpoint.addClass('clicked');

    title = corpoint.find("circle.cor").attr("title");

    // set title on contextmenu
    $("#cor-signature").text('Load ' + title);
    console.log(menuHeight);


    // Show contextmenu
    $("#cor-menu")
      .finish()
      .toggle(100) // In the right position (the mouse)
      .css({
        top: (event.pageY - menuHeight - 5) + "px",
        left: event.pageX + 5 + "px"
      });
  });
  // If the document is clicked somewhere
  $(document).on("mousedown", function(e) {
    // If the clicked element is not the menu
    if (!$(e.target).parents(".custom-menu").length > 0) {
      // Hide it
      $(".custom-menu").hide(100);
      corpoint.removeClass('clicked');

    }
  });
}

function initContextMenu(nsId) {
   // If the menu element is clicked

   $("#cor-menu li").click(function() {
    // This is the triggered action name
    switch ($(this).attr("data-action")) {
      // A case for each action. Your actions here
      case "load":
        Shiny.onInputChange(nsId, title);
        break;
      case "clear":
        alert("second");
        break;
    }

    // Hide it AFTER the action was triggered
    $(".custom-menu").hide(100);

  });
}
