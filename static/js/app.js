/**
 * Created by abarrera on 9/15/16.
 */
$(window).load(function() {
    var initial_values = {
      tp: 0,
      resolution: parseInt($('#resolution option:selected')[0].value)
    };


    // Retrieve locally stored data
    var peak_viewers = [];

    $('.gene_id').each(function(idx,elem){
        var gene_data =$(elem).data();
        var exons = utilsGGR.exonToJson(gene_data.exons),
            tps = gene_data.timepoints;

        var nextend = gene_data.nextend;
        var viewers_sel=$('.peak-viewer');
        viewers_sel.each(function(i, viewer){
            var data_mats = $(viewer).data();
            var factors = Object.keys(data_mats).filter(function(e){return e.indexOf("name")<0});
            var factor_names = Object.keys(data_mats).filter(function(e){return e.indexOf("name")>0});

            var pv = PeakviewerD3();
            var tss;
            if (exons.strand == '+')
              tss = exons.exons[0][0];
            else
              tss = exons.exons[exons.exons.length-1][1];

            pv.tss(tss);
            pv.data({
              exons: exons.exons,
              strand: exons.strand,
              tp: initial_values.tp,
              resolution: initial_values.resolution
            });

            factors.forEach(function(ee){
                var factor_name = data_mats[factor_names.filter(function(e){return ee+"name" == e})[0]];
                if (pv.data().elems)
                  pv.data().elems.push({'reads': data_mats[ee], 'name': factor_name});
                else
                  pv.data().elems = [{'reads': data_mats[ee], 'name': factor_name}];
            });
            peak_viewers.push(pv)
        });
        peak_viewers.forEach(function(e, i){
            e.container(viewers_sel[i].getAttribute("class").split(" ").map(function(e){return "."+e}).join(""));
            e.resolutions_set(gene_data.resolutions);
            e.timepoints(tps);
            e.nbins(e.data().elems[0]['reads'][0][0].length);
            e.render();
        });
        TimesliderD3.tp(tps[0]);
        TimesliderD3.timepoints(tps);
        TimesliderD3.render(".slider_container", peak_viewers);

    });

    $('#play-timeslider').on('click', function() {
        TimesliderD3.togglePlay("#play-timeslider")
    });

    $(".download-button").on('click', function() {
        var c = $($(this).siblings('.downloadable')[0])
                        .attr('class').split(" ")
                        .filter(function(e){return e.indexOf("container")>0})[0];
        downloadSVG(c);
    });

    $("#resolution").on('change', function() {
        var res = parseInt(event.target.value);
        peak_viewers.forEach(function(e){
          e.resolution(res);
          e.render();
        });
    });

});