$(document).ready(function() {
    if (window.location.search) {
        var input_params = {};
        /* process query string, e.g. ?obs=10&foo=bar */
        var params = $.map(
            window.location.search.match(/[\&\?]\w+=[^\&]+/g), 
            function(p, i) { 
                var kv = p.substring(1).split("=");
                /* NOTE: might have issue to parse some special characters here? */
                input_params[kv[0]] = decodeURIComponent(kv[1]);
            }
        );

        /* Shiny.inputBindings.getBindings() return the InputBinding instances
           for every (native) input type that Shiny supports (selectInput, textInput,
           actionButton etc.)  */
        $.each(Shiny.inputBindings.getBindings(), function(i, b) {
            /* find all inputs within a specific input type */
            var inputs = b.binding.find('#input_container');
            $.each(inputs, function(j, inp) {
                /* check if the input's id matches the key specified in the query
                   string */
                var inp_val = input_params[$(inp).attr("id")];
                if (inp_val != undefined) {
                    b.binding.setValue(inp, inp_val);
                }
            });
        });
    }
});