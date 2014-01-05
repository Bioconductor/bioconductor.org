
var bollox;


   // var loadVis = function(dataArr) {

      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {

      var arr = [];              // Handler for .ready() called.
      $.get("hits.tsv", function (data){
          bollox = data;
          var lines = data.split("\n");
          for (var i = 0;i < lines.length; i++) {
            line = lines[i];
            if (line.length ==0) continue;
            var segs = line.split("\t");
            if (i > 0) {
              segs[1] = parseInt(segs[1]);
              segs[2] = parseInt(segs[2]);
            }
            arr.push(segs);
          }


        var data = google.visualization.arrayToDataTable([
          ['Date', 'Last Year', 'This Year'],
          ['2004',  1000,      400],
          ['2005',  1170,      460],
          ['2006',  660,       1120],
          ['2007',  1030,      540]
        ]);

        var data = google.visualization.arrayToDataTable(arr);
        var options = {
          title: 'Hits for 30-day period this year and last year'
        };

        var chart = new google.visualization.LineChart(document.getElementById('chart_div'));
        chart.draw(data, options);

      });



        /*
        var data = google.visualization.arrayToDataTable([
          ['Date', 'Last Year', 'This Year'],
          ['2004',  1000,      400],
          ['2005',  1170,      460],
          ['2006',  660,       1120],
          ['2007',  1030,      540]
        ]);
*/

      }

   // }