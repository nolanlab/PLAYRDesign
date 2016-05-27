

var graphOutputBinding = new Shiny.OutputBinding();




function make_oligo_plot(svg, oligoArray, pageWidth, pageHeight, posScale)
{
    
    var barHeight = 20;
    var gap = barHeight + 4;
    
    d3.select("body").append("div").attr("id", "tag");
    makeGrid(svg, pageWidth, pageHeight, posScale);
    drawRects(svg, oligoArray, gap, barHeight, pageWidth, pageHeight, posScale);
}


function drawRects(svg, theArray, theGap, theBarHeight, w, h, posScale)
{
    
    var penalty_color_scale = d3.scale.linear()
                  .range(["red", "white"])
                  .domain(get_penalty_extent(theArray));
    
      var last_selection = null;
   
      d3.selectAll("#oligo_plot_background").on("click", function(d) {
                  if(d3.event.altKey)
                  {
                        last_selection = null;
                        Shiny.onInputChange("currently_selected_oligo", null);
                  }
            });
      
      
      
   var rectangles = svg.selectAll("oligo_line").data(theArray)
      .enter()
      .append("g")
      .attr("transform", function(d, i){
            return("translate(0," + (i*theGap) + ")");
          })
    
    
    var innerRects = rectangles.selectAll("oligo_rect")
    	.data(function(d) {return d.pair;}).enter();

    
    innerRects.append("rect")
      .attr("rx", 3)
      .attr("ry", 3)
      .attr("x", function(d){
                  return posScale(d.start);
            })
      .attr("width", function(d){
                  return (posScale(d.end) - posScale(d.start));
            })
      .attr("height", theBarHeight)
      .attr("stroke", "none")
      .attr("fill", function(d) {return(penalty_color_scale(d.penalty));});
     
    
    
    var rectText = innerRects.append("text")
      .text(function(d){
                  return d.name;
            })
      .attr("x", function(d){
                  return (((posScale(d.end) - posScale(d.start))/2) + posScale(d.start));
            })
      .attr("transform", "translate(0," + theBarHeight + ")")
      .attr("font-size", theBarHeight)
      .attr("text-anchor", "middle")
      .attr("text-height", theBarHeight)
      .attr("fill", "#000");
   
    rectText.on('mouseover', function(e) {
            var tag = "";
            var matrix = this.getScreenCTM();
            
            tag = "Name: " + e.name + "<br/>" +
                  "Unique id: " + e.unique_id + "<br/>" +
                  "Start: " + e.start + "<br/>" +
                  "End: " + e.end + "<br/>" +
                  "Penalty: " + e.penalty + "<br/>" +
                  "Sequence: " + e.sequence + "<br/>";
            var output = document.getElementById("tag");
                
                //var x = (d3.event.pageX + 8) + "px";
                //var y = (d3.event.pageY + 8) + "px";
                
            var x = (window.pageXOffset + matrix.e + parseFloat(d3.select(this).attr("x"))) + "px";
            var y = (window.pageYOffset + matrix.f + 2) + "px";
                
            output.innerHTML = tag;
            output.style.top = y;
            output.style.left = x;
            output.style.display = "block";
            x = posScale(e.start);
            var width = posScale(e.end) - posScale(e.start);
            d3.selectAll(".overlay_rect")
                  .attr("x", x)
                  .style("display", "")
                  .attr("width", width);            
            
            }).on('mouseout', function() {
                  var output = document.getElementById("tag");
                  output.style.display = "none";
                  d3.selectAll(".overlay_rect").style("display", "none");
            }).on("click", function(d) {
                  
                  if(d3.event.altKey)
                  {
                        if(!last_selection)
                        {
                              last_selection = d.unique_id;
                              Shiny.onInputChange("currently_selected_oligo", d.unique_id);
                        }
                        else
                        {
                               console.log(last_selection + "_" + d.name);
                               Shiny.onInputChange("click_select_oligo", last_selection + "_" + d.unique_id);
                               last_selection = null;
                               Shiny.onInputChange("currently_selected_oligo", null);
                        }
                  }
                  else
                        Shiny.onInputChange("click_select_oligo", d.name);
            });
}


function makeGrid(svg, w, h, posScale){
      
      svg.append("rect")
            .attr("width", w)
            .attr("height", h - 50)
            .attr("id", "oligo_plot_background")
            .attr("fill", "white")   
    
    var xAxis = d3.svg.axis()
        .scale(posScale)
        .orient('bottom')
    
    
    
    
    var grid = svg.append('g')
        .attr('class', 'x axis')
        .attr('transform', 'translate(0, ' + (h - 50) + ')')
        .call(xAxis);
}

function convert_data_for_plotting(data, x_name, y_name)
{
      var res = [];
      for (var i = 0; i < data[x_name].length; i++)
      {
            res.push({"x": data[x_name][i], "y": data[y_name][i]});
      }
      return(res);
}
      


function convert_primer3_res(data)
{
      var ret = [];      
      
      for (var i = 0; i < data.primer3.id.length; i++)
      {      
            var id = data.primer3.id[i];
            var pair_row_id = data.primer3.pair_row_id[i];
            if(!(ret[pair_row_id]))
                  ret[pair_row_id] = {pair: []}
            var start = data.primer3.start[i];
            var end = start + data.primer3.len[i] - 1;
            
            ret[pair_row_id].pair.push({name: id, start: start, end: end, type:data.primer3.type[i], unique_id:data.primer3.unique_id[i], penalty: data.primer3.penalty[i], sequence: data.primer3.sequence[i]});
      }
      return(ret);
}

function get_penalty_extent(oligoArray)
{
      var v_max = [];
      var v_min = [];
      
      
      for(var i = 0; i < oligoArray.length; i++)
      {
            v_max.push(d3.max(oligoArray[i].pair, function(d) {return d.penalty;}));
            v_min.push(d3.min(oligoArray[i].pair, function(d) {return d.penalty;}));
      }
      return([d3.min(v_min), d3.max(v_max)]);      
}


function add_overlay_rect(svg, height)
{
      svg.append("g")
            .append("rect")
            .attr("class", "overlay_rect")
            .attr("rx", 3)
            .attr("ry", 3)
            .attr("height", height)
            .style("opacity", "0.1")
            .style("display", "none");
      
      return(svg);
}

function do_plot(data, x_scale, svg, height, width, margin, svg_id, graph_class)
{
      
      
      var y_scale = d3.scale.linear().range([height, 0]);
      y_scale.domain(d3.extent(data, function(d) { return d.y; }));
                  
                  
      var x_axis = d3.svg.axis().scale(x_scale).orient("bottom");
      var y_axis = d3.svg.axis().scale(y_scale).orient("left");
                  
      var line = d3.svg.line()
            .x(function(d) { return x_scale(d.x); })
            .y(function(d) { return y_scale(d.y); });
         
         
      svg = svg.attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .attr("id", svg_id)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            
            svg.append("g")
                  .attr("class", "x axis")
                  .attr("transform", "translate(0," + height + ")")
                  .call(x_axis);
            svg.append("g")
                  .attr("class", "y axis")
                  .call(y_axis);

       
            svg.append("path")
                  .datum(data)
                  .attr("class", graph_class)
                  .attr("d", line);
                  
                  
            svg.append("g")
                  .append("rect")
                  .attr("class", "overlay_rect")
                  .attr("rx", 3)
                  .attr("ry", 3)
                  .attr("height", height)
                  .style("opacity", "0.1")
                  .style("display", "none");
            
            
            return(svg);
}



$.extend(graphOutputBinding, {
         find: function(scope) {
         return $(scope).find('.shiny-reactiveGraph-output');
         },
         renderValue: function(el, data)
         {
               var oligoArray = convert_primer3_res(data);
               var est_res = null;
               if(data.est)
                  est_res = convert_data_for_plotting(data.est, "pos", "exon_skip");
               var seq_char = convert_data_for_plotting(data.seq_char, "pos", "tm");
               var blast_res1 = convert_data_for_plotting(data.blast, "pos", "blast1");
               var blast_res2 = convert_data_for_plotting(data.blast, "pos", "blast2");
               
            
            var margin = {top: 8, right: 50, bottom: 20, left: 50}, width = 1400 - margin.left - margin.right, height = 120 - margin.top - margin.bottom;
            
      var x_scale = d3.scale.linear().range([0, width]);
      x_scale.domain(d3.extent(blast_res1, function(d) { return d.x; }));
            
            
                  

            //remove the old graph
            var svg = d3.select(el).select("svg");
            svg.remove();
         
            $(el).html("");
         
            svg = d3.select(el).append("svg");
            do_plot(blast_res1, x_scale, svg, height, width, margin, "blast1", "line1");
            svg = d3.select(el).append("svg")
            do_plot(blast_res2, x_scale, svg, height, width, margin, "blast2", "line2");
            svg = d3.select(el).append("svg")
            do_plot(seq_char, x_scale, svg, height, width, margin, "seq_char", "line3");
            if(est_res) 
            {
                  svg = d3.select(el).append("svg")
                  do_plot(est_res, x_scale, svg, height, width, margin, "est_res", "line1");
            }
            
           
            var oligo_plot_svg = d3.select(el).append("svg");
            height = height * 4;
            oligo_plot_svg = oligo_plot_svg.attr("width", width + margin.left + margin.right)
                  .attr("height", height + margin.top + margin.bottom)
                  .attr("id", "oligo_plot")
                  .append("g")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            make_oligo_plot(oligo_plot_svg, oligoArray, width, height, x_scale);

         }
});
            

Shiny.outputBindings.register(graphOutputBinding, 'reactiveGraphbinding');

