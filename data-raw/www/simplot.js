// !preview r2d3 data=c(0.3)
//
// r2d3: https://rstudio.github.io/r2d3
//

const radius = 5;
const margin = radius + 10;

var height = 50;
const rowIndex = 0;

const xScale = d3.scaleLinear();
xScale.domain([-1, 1]).range([0, 200]);

console.log(xScale)

svg.append('line')
  .attr("x1", xScale(0))
  .attr("x2", xScale(0))
  .attr("y1", -height * rowIndex) // makes match up across rows
  .attr("y2", height)
  .style("stroke", "#ddd")
  .style("shape-rendering", "crispEdges")
  .style("stroke-width", "1px")
  .style("stroke-dasharray", "3,3")


// data text
svg
  .append("g")
  .append("text")
  .attr("x", xScale(data))
  .attr("y", height)
  .attr("class", "x text")
  .attr("dy", -2)
  .style("fill", "#ddd")
  .style("font", "11px Arial, sans-serif")
  .style("text-anchor", "middle")
  .text(data);


// draw circles
svg
  .append("g")
  .append("circle")
  .attr("cx", xScale(data))
  .attr("cy", height / 2)
  .attr("r", radius)
  .style("fill", "transparent")
  .style("stroke-width", "1.1px")
  .style("stroke", "rgba(0, 0, 0, 0.75)")

