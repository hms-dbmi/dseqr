var clicks = 0;
var titles = ['only show compounds with a clinical phase', 'show all perturbations'];
function toggleClinicalTitle(el) {
 clicks += 1;
 var newTitle = clicks % 2 === 0 ? titles[1] : titles[2];
 el.title = newTitle;
}
