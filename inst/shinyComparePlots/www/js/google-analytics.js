(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
})(window,document,'script','//www.google-analytics.com/analytics.js','ga');

ga('create', 'UA-57486113-2', 'auto');
ga('send', 'pageview');

$(document).ready(function() {
  //$('input, textarea, select').change(function() {
  //  var tmp = {"db": "a", "p": "e", "items": ["b", "c", "d"]};
  //  console.log("ignore : " + JSON.stringify(tmp));
  //  ga('send', 'event', 'ignore', JSON.stringify(tmp));
  //});
  
  // For tracking on internal (NCI/DTB) group application only.
  $("#rm-predIds").change(function() { 
  	var dataset = $("#rm-dataset").val();
  	var responseId = $("#rm-responseId").val();
  	var predIds = $("#rm-predIds").val();
  	var tmp = {"d": dataset, "r": responseId, "p": predIds};
  	//console.log("rmModels : " + JSON.stringify(tmp));
    //ga('send', 'event', 'rmModels', JSON.stringify(tmp));
  });
  
	// Table search
	$('input[type="search"]').change(function() {
		var label = "tableSearch";
		var action = "";
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Dropdowns
	$('select[id]').change(function() { 
		var id = $(this).attr("id");
		
		var label = "dropdown";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Checkboxes, sliders
	$('input[id]').change(function() { 
		var id = $(this).attr("id");
		
		var label = "checkbox_slider";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Select Multi, Radio
	$('input[type]').change(function() { 
		var id = $(this).attr("type");
		
		var label = "selectizeMulti_radio";
		var action = "";
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Selectize
	$('div[data-value]').click(function() { 
		var id = $(this).attr("data-value"); 
		
		var label = "selectize";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Tabs
	$('a[data-value]').click(function() { 
		var id = $(this).attr("data-value"); 
		
		var label = "tabs";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Downloads
	$('a[id]').click(function() { 
		var id = $(this).attr("id");

		var label = "downloads";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
	
	// Button
	$('button[id]').click(function() { 
		var id = $(this).attr("id");
		
		var label = "button";
		var action = id;
		ga('send', 'event', label, action);
		
		console.log(label + " " + action);
	});
});
