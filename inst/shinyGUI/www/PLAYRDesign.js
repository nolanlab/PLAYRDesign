

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

function convert_blast_res(data)
{
      var blast_res = [];
      console.log(data.blast)
      for (var i = 0; i < data.blast.pos.length; i++)
      {
            blast_res.push({"pos": data.blast.pos[i], "blast1": data.blast.blast1[i], "blast2" : data.blast.blast2[i]});
      }
      return(blast_res);
}


function convert_est_res(data)
{
     var est_res = [];
     for (var i = 0; i < data.est.pos.length; i++)
     {
           est_res.push({"pos": data.est.pos[i], "exon_skip": data.est.exon_skip[i]});
     }
     return(est_res);
      
}

function convert_seq_char(data)
{
      var seq_char = [];
      for (var i = 0; i < data.seq_char.pos.length; i++)
      {
            seq_char.push({"pos": data.seq_char.pos[i], "gc": data.seq_char.gc[i], "tm": data.seq_char.tm[i]});
      }
      return(seq_char);
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



$.extend(graphOutputBinding, {
         find: function(scope) {
         return $(scope).find('.shiny-reactiveGraph-output');
         },
         renderValue: function(el, data)
         {
               console.log(data);
               
               var blast_res = convert_blast_res(data);
               var oligoArray = convert_primer3_res(data);
               var seq_char = convert_seq_char(data);
               var est_res = convert_est_res(data);
               
            
            var margin = {top: 50, right: 50, bottom: 50, left: 50}, width = 1400 - margin.left - margin.right, height = 300 - margin.top - margin.bottom;
            var posScale = d3.scale.linear()
                  .range([0, width]);
            var blast1_y = d3.scale.linear()
                  .range([height, 0]);
            var blast2_y = d3.scale.linear()
                  .range([height, 0]);
            var tm_y = d3.scale.linear()
                  .range([height, 0]);
            var est_y = d3.scale.linear()
                  .range([height, 0]);
                  
           
                  
            var xAxis = d3.svg.axis()
                  .scale(posScale)
                  .orient("bottom");
            var blast1_yAxis = d3.svg.axis()
                  .scale(blast1_y)
                  .orient("left");
            var blast2_yAxis = d3.svg.axis()
                  .scale(blast2_y)
                  .orient("right");
            var est_yAxis = d3.svg.axis()
                  .scale(est_y)
                  .orient("left");
                  
            var blast1_line = d3.svg.line()
                  .x(function(d) { return posScale(d.pos); })
                  .y(function(d) { return blast1_y(d.blast1); });
            var blast2_line = d3.svg.line()
                  .x(function(d) { return posScale(d.pos); })
                  .y(function(d) { return blast2_y(d.blast2); });
            var tm_line = d3.svg.line()
                  .x(function(d) { return posScale(d.pos); })
                  .y(function(d) { return tm_y(d.tm); });
            var est_line = d3.svg.line()
                  .x(function(d) { return posScale(d.pos); })
                  .y(function(d) { return est_y(d.exon_skip); });
                  
                  
                  

            //remove the old graph
            var svg = d3.select(el).select("svg");
                  svg.remove();
         
            $(el).html("");
         
            //append a new one
            svg = d3.select(el).append("svg");
            svg = svg.attr("width", width + margin.left + margin.right)
                  .attr("height", height + margin.top + margin.bottom)
                  .attr("id", "main_graph")
                  .append("g")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            posScale.domain(d3.extent(blast_res, function(d) { return d.pos; }));
            blast1_y.domain(d3.extent(blast_res, function(d) { return d.blast1; }));
            blast2_y.domain(d3.extent(blast_res, function(d) { return d.blast2; }));
            tm_y.domain(d3.extent(seq_char, function(d) { return d.tm; }));
            
            svg.append("g")
                  .attr("class", "x axis")
                  .attr("transform", "translate(0," + height + ")")
                  .call(xAxis);
            svg.append("g")
                  .attr("class", "y axis")
                  .call(blast1_yAxis)
            svg.append("g")
                  .attr("class", "y axis")
                  .attr("transform", "translate(" + width + " ,0)")  
                  .call(blast2_yAxis)

       
            svg.append("path")
                  .datum(blast_res)
                  .attr("class", "line1")
                  .attr("d", blast1_line);
            
            svg.append("path")
                  .datum(blast_res)
                  .attr("class", "line2")
                  .attr("d", blast2_line);
                  
            svg.append("path")
                  .datum(seq_char)
                  .attr("class", "line3")
                  .attr("d", tm_line);
            
            svg = add_overlay_rect(svg, height);
            var est_plot_svg = d3.select(el).append("svg");
            
            est_plot_svg = est_plot_svg.attr("width", width + margin.left + margin.right)
                  .attr("height", height + margin.top + margin.bottom)
                  .attr("id", "est_graph")
                  .append("g")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            posScale.domain(d3.extent(est_res, function(d) { return d.pos; }));
            est_y.domain(d3.extent(est_res, function(d) { return d.exon_skip; }));
            
            est_plot_svg.append("g")
                  .attr("class", "x axis")
                  .attr("transform", "translate(0," + height + ")")
                  .call(xAxis);
            est_plot_svg.append("g")
                  .attr("class", "y axis")
                  .call(est_yAxis);
       
            est_plot_svg.append("path")
                  .datum(est_res)
                  .attr("class", "line1")
                  .attr("d", est_line);
            est_plot_svg = add_overlay_rect(est_plot_svg, height);
            var oligo_plot_svg = d3.select(el).append("svg");
            height = height * 2;
            oligo_plot_svg = oligo_plot_svg.attr("width", width + margin.left + margin.right)
                  .attr("height", height + margin.top + margin.bottom)
                  .attr("id", "oligo_plot")
                  .append("g")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            make_oligo_plot(oligo_plot_svg, oligoArray, width, height, posScale);


         }
});
            

Shiny.outputBindings.register(graphOutputBinding, 'reactiveGraphbinding');

