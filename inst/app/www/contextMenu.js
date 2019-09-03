var title = null;

function setupContextMenu() {
  title = null;

  $('.cor').off().on("contextmenu", function(event) {
    // Avoid the real one
    event.preventDefault();

    title = $(event.target).attr("title");

    // Show contextmenu
    $("#cor-menu")
      .finish()
      .toggle(100) // In the right position (the mouse)
      .css({
        top: event.pageY + "px",
        left: event.pageX + "px"
      });
  });
  // If the document is clicked somewhere
  $(document).off().on("mousedown", function(e) {
    // If the clicked element is not the menu
    if (!$(e.target).parents(".custom-menu").length > 0) {
      // Hide it
      $(".custom-menu").hide(100);
    }
  });
}

function initContextMenu() {
   // If the menu element is clicked
   $("#cor-menu li").click(function() {
    // This is the triggered action name
    switch ($(this).attr("data-action")) {
      // A case for each action. Your actions here
      case "load":
        console.log(title);
        break;
      case "clear":
        alert("second");
        break;
    }

    // Hide it AFTER the action was triggered
    $(".custom-menu").hide(100);
  });
}
