Shiny.addCustomMessageHandler("showLoading",
	function(message) {
		//console.log(message);
		
		if(message.show) {
			$("#loadingModal").modal("show");
		} else {
			$("#loadingModal").modal("hide");
		}
	}
);