import os
import pysam
import numpy




def create_header(bin_list, directory, active, long_read_qc_html, short_read_qc_html=None):

    main_header = '''<!doctype html>
<html lang="en">
	<head>
        <title>SDMass</title>

        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">


        <!-- Main CSS -->
        <link rel="stylesheet" href="''' + directory + '''css/style.css">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.css">
        <!-- Font Awesome -->
        <link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet">
        <script
			  src="https://code.jquery.com/jquery-3.3.1.min.js"
			  integrity="sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8="
			  crossorigin="anonymous"></script>
	    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
	    <script type="text/javascript" class="init">

$(document).ready(function() {
	$('#thetable').DataTable();
} );

        </script>
    </head>


    <body>
  
        
        <!-- Main navigation -->
        <div id="sidebar">
                      
            <div class="navbar-expand-md navbar-dark"> 
            
                <header class="d-none d-md-block">
                    <h1><span>SD</span>Mass</h1>
                </header>
                
                
                <!-- Mobile menu toggle and header -->
                <div class="mobile-header-controls">
                    <a class="navbar-brand d-md-none d-lg-none d-xl-none" href="#"><span>SD</span>Mass</a>
                    <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#SidebarContent" aria-controls="SidebarContent" aria-expanded="false" aria-label="Toggle navigation">
                        <span class="navbar-toggler-icon"></span>
                    </button>
                </div>
         
                <div id="SidebarContent" class="collapse flex-column navbar-collapse">
 
                        
                    
                    <!-- Main navigation items -->
                    <nav class="navbar navbar-dark">
                        <div id="mainNavbar">
                            <ul class="flex-column mr-auto">
'''
    if active == 'index.html':
        main_header += '                                <li class="nav-item active">\n'
    else:
        main_header += '                                <li class="nav-item">\n'

    main_header += '''                                    <a class="nav-link" href="''' + directory + '''index.html">Home <span class="sr-only">(current)</span></a>
                                </li>

                                <li class="nav-item">
                                        <a class="nav-link" href="''' + directory + long_read_qc_html + '''">Long read stats</a>
                                </li>'''

    if not short_read_qc_html is None:
        main_header += '''                                <li class="nav-item">
                                        <a class="nav-link" href="''' + directory + short_read_qc_html + '''">Short read stats</a>
                                </li>
'''
    main_header += '''                                <li class="nav-item dropdown">
                                            <a class="nav-link dropdown-toggle" href="#MenuDropdown" data-toggle="collapse" aria-controls="MenuDropdown" aria-expanded="false">Details of bins</a>
                                            <ul id="MenuDropdown" class="sub-navbar collapse flex-column">'''
    for i in bin_list:
        if active == i:
            main_header += '                                                <li class="nav-item active"><a class="nav-link" href="%sbin/%s.html">Bin %s</a></li>\n' % (directory, i, i)
        else:
            main_header += '                                                <li class="nav-item"><a class="nav-link" href="%sbin/%s.html">Bin %s</a></li>\n' % (directory, i, i)
    main_header += '''                                            </ul>
                                </li>
                                <li class="nav-item">
                                        <a class="nav-link" href="gtdbtk.html">GTDBtk</a>
                                </li>
                                <li class="nav-item">
                                        <a class="nav-link" href="https://github.com/mjsull/SDMass/issues">Help</a>
                                </li>
                            </ul>
                        </div>   
                    </nav>
                
                </div>
            </div> 
        </div>        
        
        
        <div id="content">
            <div id="content-wrapper">
'''


    return(main_header)

def add_title(header, subheader=''):
    return('''
            <!-- Jumbtron / Slider -->
                <div class="jumbotron-wrap">
                    <div class="container-fluid">
                        <div class="jumbotron static-slider">
                            <h1 class="text-center">''' + header + \
    '''</h1>
                            <p class="lead text-center">''' + subheader + '''</p>
                        </div>
                    </div>
                </div>
''')

def add_main(header, text, contigs=None):
    if contigs is None:
        return('''                <main class="container-fluid">
                    <div class="row">
                        <!-- Main content -->
                        <div class="col-md-12">
                            <article>
                                <h2 class="article-title">''' + header + '''</h2>

                                <p> ''' + text + ''' </p>
''')
    else:
        html_string = '''                <main class="container-fluid">
                    <div class="row">
                        <!-- Sidebar -->
                        <aside class="col-md-2">
                            <div class="sidebar-box">
                                <h4>Contigs</h4>
                                <div class="list-group list-group-root">
                                    <a class="list-group-item active" href="index.html">Overview</a>
'''
        count = 0
        for i in contigs:
            the_path = "www/contigs/" + i + ".html"
            if os.path.exists(the_path):
                html_string += '                                    <a class="list-group-item" href="' + the_path + '">'+ i + '</a>\n'
            else:
                count += 1
        html_string += '''                                </div>
                            </div>

                            <div class="sidebar-box sidebar-box-bg">
                                <h4>n.b.</h4>
                                <p> ''' + str(count) + ''' contigs under 100Kbp not shown.</p>
                            </div>

                        </aside>

                        <!-- Main content -->
                        <div class="col-md-10">
                            <article>
                                <h2 class="article-title">''' + header + '''</h2>

                                <p> ''' + text + ''' </p>
'''
        return(html_string)

def end_main():
    return('''


                        </div>





                    </div> 
                </main>
''')

def add_footer():
    return('''                <!-- Footer -->
                <div class="container-fluid footer-container">
                    <footer class="footer">
                        <div class="footer-bottom">
                                <p class="text-center">Created by <a href="https://mjsull.github.io"</a>mjsull.github.io</p>
                                <p class="text-center"><a href="#">Back to top</a></p>
                        </div>
                    </footer>
                </div> 
            </div>
        </div>

        <!-- Bootcamp JavaScript -->
        <!-- jQuery first, then Popper.js, then Bootstrap JS -->


    </body>
</html>''')

def create_table(headers, list_of_vals):
    out_string = '''                                <table class="table" id="thetable" style="width:100%">
                                    <thead>
                                        <tr>
'''
    for i in headers:
        out_string += '                                            <th>' + i + '</th>\n'
    out_string += '''                                        </tr>
                                    </thead>
                                    <tbody>
'''
    for i in list_of_vals:
        out_string += '                                        <tr>\n'
        for j in i:
            out_string += '                                            <td>' + str(j) + '</td>\n'
        out_string += '                                        </tr>\n'
    out_string += '''                                    </tbody>
                                </table>
'''
    return(out_string)





def get_cov_stats_long(bamfile, contig, bin_size=3000, bin_step=500, buffer=50):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_length = samfile.get_reference_length(contig)
    bin_num = ref_length // bin_step +1
    coverage_forward = numpy.zeros(bin_num)
    coverage_reverse = numpy.zeros(bin_num)
    trimmed_starts = numpy.zeros(bin_num, dtype=int)
    trimmed_ends = numpy.zeros(bin_num, dtype=int)
    through = numpy.zeros(bin_num, dtype=int)
    starts_in = numpy.zeros(bin_num, dtype=int)
    ends_in = numpy.zeros(bin_num, dtype=int)
    x = [y for y in range(0, ref_length+1, bin_num)]
    for read in samfile.fetch(contig):
        rstart = read.reference_start
        rend = read.reference_end
        qstart = read.query_alignment_start
        qend = read.query_alignment_end
        qlength = read.infer_query_length()
        if qstart > buffer:
            trimmed_start = True
        else:
            trimmed_start = False
        if qend < qlength - buffer:
            trimmed_end = True
        else:
            trimmed_end = False
        if read.is_secondary:
            continue
        for i in range(max([0, (rstart-bin_size) //bin_step]), rend//bin_step+1):
            if i * bin_step <= rstart < rend < i*bin_step + bin_size:
                pass
            elif i * bin_step <= rstart < i*bin_step + bin_size:
                if trimmed_start:
                    trimmed_starts[i] += 1
                else:
                    starts_in[i] += 1
            elif i * bin_step <= rend < i * bin_step + bin_size:
                if trimmed_end:
                    trimmed_ends[i] += 1
                else:
                    ends_in[i] += 1
            elif rstart <= i * bin_step and rend >= i * bin_step + bin_size:
                through[i] += 1
            bases = (min([rend, i*bin_step + bin_size]) - max([rstart, i*bin_step])) / bin_size
            if read.is_reverse:
                coverage_reverse[i] += bases
            else:
                coverage_forward[i] += bases
    return(coverage_forward, coverage_reverse, trimmed_starts, trimmed_ends, starts_in, ends_in, x)


def get_gtdbtk(gtdbtk_folder, in_dict=None):
    connect_dict = {}
    out_dict = {}
    with open(os.path.join(gtdbtk_folder, 'gtdbtk.bac120.summary.tsv')) as f:
        f.readline()
        for line in f:
            the_bin, phylo, nearest, ani_radius, ani_tax, ani = line.split('\t')[:6]
            the_bin = the_bin.split('.')[-1]
            out_dict[the_bin] = (nearest, ani)
            if not in_dict is None:
                cov = in_dict[the_bin]
            else:
                cov = 10
            lastname = None
            for i in phylo.split(';'):
                if i == 's__':
                    the_name = 's__' +the_bin
                elif i == 'g__':
                    the_name = 'g__' + the_bin
                else:
                    the_name = i
                if not lastname is None:
                    if (lastname, the_name) in connect_dict:
                        connect_dict[(lastname, the_name)] += cov
                    else:
                        connect_dict[(lastname, the_name)] = cov
                lastname = the_name
    with open(os.path.join(gtdbtk_folder, 'gtdbtk.ar122.summary.tsv')) as f:
        f.readline()
        for line in f:
            the_bin, phylo, nearest, ani_radius, ani_tax, ani = line.split('\t')[:6]
            out_dict[the_bin] = (nearest, ani)
            if not in_dict is None:
                cov = in_dict[the_bin]
            else:
                cov = 10
            lastname = None
            for i in phylo.split(';'):
                if i == 's__':
                    the_name = 's__' + the_bin
                elif i == 'g__':
                    the_name = 'g__' + the_bin
                else:
                    the_name = i
                if not lastname is None:
                    if (lastname, the_name) in connect_dict:
                        connect_dict[(lastname, the_name)] += cov
                    else:
                        connect_dict[(lastname, the_name)] = cov
                lastname = the_name
    if in_dict is None:
        return out_dict
    connect_list = []
    for i in connect_dict:
        connect_list.append([i[0], i[1], connect_dict[i]])
    connect_list.sort(key=lambda x: x[2], reverse=True)
    with open('www/gtdbtk.html', 'w') as o:
        o.write('''<html>
<body>
 <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>

<div id="sankey_multiple" style="width: 800px; height: 1200px;"></div>

<script type="text/javascript">
  google.charts.load("current", {packages:["sankey"]});
  google.charts.setOnLoadCallback(drawChart);
   function drawChart() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'From');
    data.addColumn('string', 'To');
    data.addColumn('number', 'Weight');
    data.addRows([\n''')
        for i in connect_list[:-1]:
            o.write(str(i) + ',\n')
        o.write(str(connect_list[-1]))
        o.write('''
    ]);

    // Set chart options
    var options = {
      width: 1400,
      height: 2000, 
      textStyle: {
            fontSize: 6,
        },
    };

    // Instantiate and draw our chart, passing in some options.
    var chart = new google.visualization.Sankey(document.getElementById('sankey_multiple'));
    chart.draw(data, options);
   }
</script>
</body>
</html>''')


def get_cov_stats_short(bamfile, contig, bin_size=3000, bin_step=500):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_length = samfile.get_reference_length(contig)
    bin_num = ref_length // bin_step +1
    coverage_forward = numpy.zeros(bin_num)
    coverage_reverse = numpy.zeros(bin_num)
    for read in samfile.fetch(contig):
        rstart = read.reference_start
        rend = read.reference_end
        if read.is_secondary or rstart is None or rend is None:
            continue
        for i in range(max([0, (rstart-bin_size)//bin_step]), rstart//bin_step+1):
            bases = (min([rend, i*bin_step + bin_size]) - max([rstart, i*bin_step])) / bin_size
            if read.is_reverse:
                coverage_reverse[i] += bases
            else:
                coverage_forward[i] += bases
    return coverage_forward, coverage_reverse


def get_gene_sizes(gff_file):
    size_dict = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            contig, prod, feat, start, stop = line.split()[:5]
            if not contig in size_dict:
                size_dict[contig] = []
            size_dict[contig].append(int(stop) - int(start))
    return size_dict






# def create_contig_page(bin, ctg, cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in, ends_ind,
#                        cov_forward_ill, cov_reverse_ill, x, outpage):
#      html_out = open(out_directory + '/qc_website/ctgs/' + str(html_name) + '_graphs.html', 'w')
#      html_out.write(header)
#      html_out.write('  <script type="text/javascript">\n'
#                     '  window.onload = function () {\n'
#                     '  var chart1 = new CanvasJS.Chart("chartContainer1",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Total number of each read type in bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16,\n')
# #         if clipped_flag != [] or unclipped_flag != []:
# #             html_out.write('      stripLines:[\n')
# #             for values in clipped_flag[:-1]:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(values[0]) + ',\n'
# #                                                                          '        endValue: ' + str(
# #                     values[1]) + ',\n'
# #                                  '        color: "#16CBEF"\n'
# #                                  '      },\n')
# #             if clipped_flag != []:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(clipped_flag[-1][0]) + ',\n'
# #                                                                                    '        endValue: ' + str(
# #                     clipped_flag[-1][1]) + ',\n'
# #                                            '        color: "#16CBEF"\n'
# #                                            '      }')
# #                 if unclipped_flag != []:
# #                     html_out.write(',')
# #                 else:
# #                     html_out.write('\n')
# #             for values in unclipped_flag[:-1]:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(values[0]) + ',\n'
# #                                                                          '        endValue: ' + str(
# #                     values[1]) + ',\n'
# #                                  '        color: "#53A01F"\n'
# #                                  '      },\n')
# #             if unclipped_flag != []:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(unclipped_flag[-1][0]) + ',\n'
# #                                                                                      '        endValue: ' + str(
# #                     unclipped_flag[-1][1]) + ',\n'
# #                                              '        color: "#53A01F"\n'
# #                                              '      }\n')
# #             html_out.write('       ]\n  ')
# #
#     html_out.write('      },\n'
#                    '      axisY: {\n'
#                    '      title: "Number of reads",\n'
#                    '      titleFontSize: 16\n'
#                    '      },\n'
#                    '       data: [')
#     for the_data, leg_lab, colour in zip(
#             [forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end,
#              forward_start_clipped,
#              reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
#                                                                                 'Read spans bin (Reverse)',
#                                                                                 'Read starts in bin (F)',
#                                                                                 'Read starts in bin (R)',
#                                                                                 'Read terminates in bin (F)',
#                                                                                 'Read terminates in bin (R)',
#                                                                                 'Read start clipped in bin (F)',
#                                                                                 'Read start clipped in bin (R)',
#                                                                                 'Read end clipped in bin (F)',
#                                                                                 'Read end clipped in bin (R)'],
#             ["#BE4200", "#EE16B8", "#38081F", "#4282F6", "#45FAA5", "#763181",
#              "#234D65", "#F1F272", "#D15AF4", "#CBC22F"]):
#         html_out.write('\n{\n'
#                        '         type: "line",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         color: "' + colour + '",\n'
#                                                                                            '         markerType: "none",\n'
#                                                                                            '         dataPoints: [\n')
#         if refnum == 0:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'w')
#             bgt_out = open(
#                 options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.bwt',
#                 'w')
#             bgt_out.write(
#                 'track type=bigwig bigDataUrl=https://vanbah01.u.hpc.mssm.edu/igb/' + options.assembly_name + '/'
#                 + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower()
#                 + '.bw name=' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() +
#                 'color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
#             bgt_out.close()
#         else:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'a')
#         bg_out.write(
#             'fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#             if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
#                 bg_out.write(str(the_data[value]) + '\n')
#         bg_out.close()
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is reverse_end_clipped:
#             html_out.write(',\n')
#     html_out.write('      ],\n'
#                    '      rangeChanged: syncHandler\n'
#                    '    });\n'
#                    '    chart1.render();\n'
#                    '  var chart2 = new CanvasJS.Chart("chartContainer2",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Proportion of each read type in bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16\n'
#                                                      '},\n'
#                                                      '      axisY: {\n'
#                                                      '      title: "Number of reads",\n'
#                                                      '      titleFontSize: 16\n'
#                                                      '},\n'
#                                                      '       data: [')
#     for the_data, leg_lab in zip(
#             [forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end,
#              forward_start_clipped,
#              reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
#                                                                                 'Read spans bin (Reverse)',
#                                                                                 'Read starts in bin (F)',
#                                                                                 'Read starts in bin (R)',
#                                                                                 'Read terminates in bin (F)',
#                                                                                 'Read terminates in bin (R)',
#                                                                                 'Read start clipped in bin (F)',
#                                                                                 'Read start* clipped in bin (R)',
#                                                                                 'Read end clipped in bin (F)',
#                                                                                 'Read end clipped in bin (R)']):
#         html_out.write('         {\n'
#                        '         type: "stackedArea100",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         markerType: "none",\n'
#                                                             '         legendMarkerType: "square",\n'
#                                                             '        dataPoints: [\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is reverse_end_clipped:
#             html_out.write(',')
#     html_out.write('      ],\n'
#                    '      rangeChanged: syncHandler\n'
#                    '   });\n'
#                    '    chart2.render();\n'
#                    '  var chart3 = new CanvasJS.Chart("chartContainer3",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Number of large indels in each bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + '",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16,\n')
#     if indel_flag != []:
#         html_out.write('      stripLines:[\n')
#         for values in indel_flag[:-1]:
#             html_out.write('      {\n'
#                            '        startValue: ' + str(values[0]) + ',\n'
#                                                                      '        endValue: ' + str(
#                 values[1]) + ',\n'
#                              '        color: "#16CBEF"\n'
#                              '      },\n')
#         html_out.write('      {\n'
#                        '        startValue: ' + str(indel_flag[-1][0]) + ',\n'
#                                                                          '        endValue: ' + str(
#             indel_flag[-1][1]) + ',\n'
#                                  '        color: "#16CBEF"\n'
#                                  '      }\n'
#                                  '      ]\n')
#     html_out.write('\n    },\n'
#                    '      axisY: {\n'
#                    '      title: "Number of reads",\n'
#                    '      titleFontSize: 16\n'
#                    '},\n'
#                    '      data: [')
#     for the_data, leg_lab in zip([large_deletions, large_insertions, coverage_array],
#                                  ['Deletions in read', 'Insertions in read', 'Total reads']):
#         html_out.write('         {\n'
#                        '         type: "line",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         markerType: "none",\n'
#                                                             '        dataPoints: [\n')
#         if refnum == 0:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'w')
#             bgt_out = open(
#                 options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.bwt',
#                 'w')
#             bgt_out.write(
#                 'track type=bigwig bigDataUrl=' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                      '').lower()
#                 + '.bw name=test color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
#             bgt_out.close()
#         else:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'a')
#         bg_out.write(
#             'fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#             if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
#                 bg_out.write(str(the_data[value]) + '\n')
#         bg_out.close()
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is coverage_array:
#             html_out.write(',')
#     html_out.write('''      ],
#     rangeChanged: syncHandler\n
# });
#
# chart3.render();
# var charts = [chart1, chart2, chart3];
#
# function syncHandler(e) {
#
# for (var i = 0; i < charts.length; i++) {
#     var chart = charts[i];
#
#     if (!chart.options.axisX) chart.options.axisX = {};
#
#     if (e.trigger === "reset") {
#         chart.options.axisX.viewportMinimum = chart.options.axisX.viewportMaximum = null;
#
#     } else if (chart !== e.chart) {
#         chart.options.axisX.viewportMinimum = e.axisX.viewportMinimum;
#         chart.options.axisX.viewportMaximum = e.axisX.viewportMaximum;
#     }
#
#     chart.render();
#
# }
# }
# function clickHandler(e) {
# var x = parseInt(e.target.id.split(',')[0])
# var y = parseInt(e.target.id.split(',')[1])
# for (var i = 0; i < charts.length; i++) {
#     var chart = charts[i];
#     chart.options.axisX.viewportMinimum = x;
#     chart.options.axisX.viewportMaximum = y;
#     chart.render();
#     }
# }
# var zoomButtons = document.getElementsByClassName("zoom");
# for (var i = 0; i < zoomButtons.length; i++) {
# var zoomButton = zoomButtons[i]
# zoomButton.addEventListener("click", clickHandler)
# };
# }
# </script>
# <script type="text/javascript" src="/igb/webpage_css_js/canvasjs.min.js"></script></head>
# <h1> ''' + str(html_name) + ''' graphs </h1>
# <div id="chartContainer1" style="height: 400px; width: 100%;">
# </div>
# <br />
# <table style="width: 60%" align="center">
#     <tr>
#         <th>Type</th>
#         <th>Position</th>
#     </tr>
# ''')
#     for j in clipped_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     for j in unclipped_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     if unclipped_flag == [] and clipped_flag == []:
#         html_out.write('        <tr>\n          <td> no flags </td>\n')t.write('          <td> no flags </td>\n        </tr>\n')
#     html_out.write('''
# </table>
# <br />
# <div id="chartContainer2" style="height: 400px; width: 100%;">
# </div>
# <br />
# <div id="chartContainer3" style="height: 400px; width: 100%;">
# </div>
# <br />
# <table style="width: 60%" align="center">
#     <tr>
#         <th>Type</th>
#         <th>Position</th>
#     </tr>
# ''')
#     for j in indel_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     if indel_flag == []:
#         html_out.write('        <tr>\n          <td> no flags </td>\n')
#         html_out.write('          <td> no flags </td>\n        </tr>\n')
#     html_out.write('''
# </table>
# <br />
# ''')
#     html_out.write(footer)
#     html_out.close()
#     if sam.lengths[refnum] >= 300000:
#         indexa = 100000 / bin_step
#     else:
#         indexa = sam.lengths[refnum] / 3 / bin_step
#     new_array = coverage_array[indexa:-indexa]
#     if len(new_array) >= 1000:
#         zesteps = len(new_array) / 1000
#     else:
#         zesteps = 1
#     covlist = []
#     for jump in range(0, len(new_array) + 1, zesteps):
#         covlist.append(sum(new_array[jump:jump + zesteps]) * 1.0 / zesteps)
#     out_cov[html_name] = covlist
#     out_flag[html_name] = (len(clipped_flag), len(unclipped_flag), len(indel_flag))
#     csfile = open(options.output_folder + '/chrom.size', 'w')
#     for i in chrom_size:
#         csfile.write(i + '\t' + str(chrom_size[i]) + '\n')
#     csfile.close()
#     for i in os.listdir(options.output_folder + '/wiggle/'):
#         subprocess.Popen(
#             'wigToBigWig ' + options.output_folder + '/wiggle/' + i + ' ' + options.output_folder + '/chrom.size ' + options.output_folder + '/bigwig/' + i[
#                                                                                                                                                           :-3] + 'bw',
#             shell=True).wait()
#     return out_cov, out_flag


def create_main_page(outfile, fasta, checkm_file, metabat_folder, long_bam, short_bam, gff_file, long_qc_html, short_qc_html, gtdbtk_dir, min_contig_size=100000):
    with open(fasta) as f:
        len_dict = {}
        for line in f:
            if line.startswith('>'):
                name = line.split()[0][1:]
                len_dict[name] = 0
            else:
                len_dict[name] += len(line.rstrip())
    checkm_dict = {}
    with open(checkm_file) as f:
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            if line.startswith('-----'):
                break
            bin_id, marker_lineage, genomes, markers, marker_sets, e0, e1, e2, e3, e4, e5, completeness, contamination, heterogeneity = [s.strip() for s in line.rstrip().split('  ') if s]
            checkm_dict[bin_id.split('.')[1]] = [marker_lineage, completeness, contamination, heterogeneity]
    outlist = []
    bin_list = []
    gene_size_dict = get_gene_sizes(gff_file)
    for i in os.listdir(metabat_folder):
        if not i.startswith('binned_contigs'):
            continue
        bin = i.split('.')[1]
        bin_list.append(bin)
    gtdbtk_dict = get_gtdbtk(gtdbtk_dir)
    for i in gtdbtk_dict:
        gtdbtk_dict[i] = gtdbtk_dict[i][0] + ' (' + gtdbtk_dict[i][1] + '%)'
    cov_dict = {}
    for i in os.listdir(metabat_folder):
        if not i.startswith('binned_contigs'):
            continue
        bin = i.split('.')[1]
        ctgs = []
        with open(os.path.join(metabat_folder, i)) as fa:
            for line in fa:
                if line.startswith('>'):
                    ctg_name = line.split()[0][1:]
                    ctgs.append(ctg_name)
        length_list = []
        for j in ctgs:
            length_list.append(len_dict[j])
        max_contig = max(length_list)
        bases_assembled = sum(length_list)

        length_list.sort(reverse=True)
        y = 0
        for j in length_list:
            y += j
            if y >= bases_assembled/2:
                break
        n50 = j
        ctg_details = []
        bases_sequenced_long = 0
        bases_sequenced_short = 0
        gene_sizes = []
        for ctg in ctgs:
            gene_sizes += gene_size_dict[ctg]
            cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in, ends_ind, x = get_cov_stats_long(long_bam, ctg)
            coverage_long = (sum(cov_forward) + sum(cov_reverse)) / len(cov_forward)
            bases_sequenced_long += coverage_long * len_dict[ctg]
            if short_bam is None:
                cov_forward_ill, cov_reverse_ill = None, None
                coverage_short = 0
            else:
                cov_forward_ill, cov_reverse_ill = get_cov_stats_short(short_bam, ctg)
                coverage_short = (sum(cov_forward_ill) + sum(cov_reverse_ill)) / len(cov_forward_ill)
                bases_sequenced_short += coverage_short * len_dict[ctg]
            # if len_dict[ctg] >= min_contig_size:
            #     create_contig_page(bin, ctg, cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in,
            #                     ends_ind, cov_forward_ill, cov_reverse_ill, x, 'ctg/' + ctg + '.html')
            ctg_details.append((ctg, len_dict[ctg], coverage_long, coverage_short))
        gene_average = numpy.average(gene_sizes)
        gene_std = numpy.std(gene_sizes)
        gene_no = len(gene_sizes)
        if short_bam is None:
            cov_dict[bin] = bases_sequenced_long/bases_assembled
        else:
            cov_dict[bin] = bases_sequenced_short/bases_assembled
        bin_headers = ['Bin', 'Max. contig (bp)', '# of contigs', 'bases assembled', 'N50', 'average read depth (long)',
                              'average read depth (short)', 'average gene size', 'Gene size Std. dev.', '# of genes', 'marker lineage',
                              'completeness', 'Contamination', 'Heterozygosity', 'Best mash hit']
        bin_details = [bin, '{:,}'.format(max_contig), '{:,}'.format(len(ctgs)), '{:,}'.format(bases_assembled), '{:,}'.format(n50),
                        '{:,.2f}'.format(bases_sequenced_long/bases_assembled),'{:,.2f}'.format(bases_sequenced_short/bases_assembled),
                       '{:,.2f}'.format(gene_average), '{:,.2f}'.format(gene_std), '{:,}'.format(gene_no),
                       ] + checkm_dict[bin] + [gtdbtk_dict[bin]]
        create_bin_page(bin_headers, bin_details, ctg_details, 'bin/' + bin + '.html', bin_list, long_qc_html, short_qc_html)
        outlist.append(bin_details)
    get_gtdbtk(gtdbtk_dir, cov_dict)
    with open(outfile, 'w') as o:
        o.write(create_header(bin_list, '', 'index.html', long_qc_html, short_qc_html))
        o.write(add_title("Overview of assembly", "assembled into " + str(len(bin_list)) + " bins."))
        o.write(add_main("Overview of bins", "details for each bin"))
        o.write(create_table(bin_headers, outlist))
        o.write(end_main())
        o.write(add_footer())







    # write page # of bins
    # for each bin
         # get max contig, bases assembled, N50,
         # get average read depth (long)
         # get average read depth (short)
         # get marker lineage, completeness, Contamination, Heterozygosity
         # get best mash hit

         # add entry to table (Bin 	Max. contig 	bases assembled 	N50 	average read depth (long) 	average read depth (short) 	marker lineage 	complete 	Contamination 	Heterozygosity 	Best mash hit)
         # generate bin page





def create_bin_page(headers, bin_details, ctg_details, outfile, bin_list, long_qc_html, short_qc_html):
    headers = ['Bin', 'Max. contig (bp)', '# of contigs', 'bases assembled', 'N50', 'average read depth (long)',
                              'average read depth (short)', 'marker lineage',  'completeness', 'Contamination',
                              'Heterozygosity', 'Best mash hit']
    with open('www/' + outfile, 'w') as o:
        o.write(create_header(bin_list, '../', bin_details[0], long_qc_html, short_qc_html))
        o.write(add_title("overview of bin " + bin_details[0], "assembled into " + bin_details[2] + " contigs."))
        main_string = ''
        for i, j in zip(headers, bin_details):
            main_string +=  '<b>' + i + ':</b> ' + j + '<br>\n'
        ctgs = []
        for i in ctg_details:
            ctgs.append(i[0])
        o.write(add_main("Bin details:", main_string, ctgs))
        o.write(create_table(["contig", "length", "coverage long", "coverage short"], ctg_details))
        o.write(end_main())
        o.write(add_footer())



    # create table with contig/coverage long/ coverage short/length/mash hit?




checkm_file = snakemake.input.checkm_file
metabat_bins = snakemake.input.metabat_bins[:-5]
fasta = snakemake.input.fasta
long_html = snakemake.input.long_reads_qc_html[4:]
short_html = snakemake.input.short_reads_qc_html[4:]
gff_file = snakemake.input.genes_gff
gtdbtk_dir = snakemake.input.gtdbtk_summary[:-25]
try:
    os.makedirs('www/bin')
except FileExistsError:
    pass



if os.path.exists("data/short_reads.fastq.gz"):
    short_bam = "data/final_short.sort.bam"
    long_bam = "data/final_long.sort.bam"
else:
    long_bam = "data/final_long.sort.bam"
    short_bam = None


create_main_page("www/index.html", fasta, checkm_file, metabat_bins, long_bam, short_bam, gff_file, long_html, short_html, gtdbtk_dir)