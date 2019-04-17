import os
import pysam
import numpy




def create_header(bin_list, active):

    main_header = '''<!doctype html>
<html lang="en">
	<head>
        <title>SDMass</title>

        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">


        <!-- Main CSS -->
        <link rel="stylesheet" href="css/style.css">
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
                                <li class="nav-item">
                                        <a class="nav-link" href="index.html">Home <span class="sr-only">(current)</span></a>
                                </li>

                                <li class="nav-item">
                                        <a class="nav-link" href="long_reads.html">Long read stats</a>
                                </li>

                                <li class="nav-item">
                                        <a class="nav-link" href="short_reads.html">Short read stats</a>
                                </li>

                                <li class="nav-item dropdown">
                                            <a class="nav-link dropdown-toggle" href="#MenuDropdown" data-toggle="collapse" aria-controls="MenuDropdown" aria-expanded="false">Details of bins</a>
                                            <ul id="MenuDropdown" class="sub-navbar collapse flex-column">'''
    for i in bin_list:
        main_header += '                                                <li class="nav-item"><a class="nav-link" href="bin%s.html">Bin %s</a></li>\n' (bin, bin)
    main_header += '''                                            </ul>
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
    main_header.replace('<li class="nav-item">\n                                        <a class="nav-link" href="' + active +
                        '">Home <span class="sr-only">(current)</span></a>',
                        '<li class="nav-item active">\n                                        <a class="nav-link" href="'
                        + active + '">Home <span class="sr-only">(current)</span></a>')

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

def add_main(header, text, width):
    return('''                <main class="container-fluid">
                    <div class="row">
                        <!-- Main content -->
                        <div class="col-md-''' + width + ''''">
                            <article>
                                <h2 class="article-title">''' + header + ''''</h2>

                                <p> ''' + text + ''' </p>
''')

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
            out_string += '                                            <td>' + str(i) + '</td>\n'
        out_string += '                                        </tr>\n'
    out_string += '''                                    </tbody>
                                </table>
'''
    return(out_string)





def get_cov_stats_long(bamfile, contig, bin_size=3000, bin_step=500, buffer=50):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_length = samfile.get_reference_length(contig)
    bin_num = ref_length / bin_step +1
    coverage_forward = numpy.zeros(bin_num, dtype=float)
    coverage_reverse = numpy.zeros(bin_num, dtype=float)
    trimmed_starts = numpy.zeros(bin_num, type=int)
    trimmed_ends = numpy.zeros(bin_num, type=int)
    through = numpy.zeros(bin_num, type=int)
    starts_in = numpy.zeros(bin_num, type=int)
    ends_in = numpy.zeros(bin_num, type=int)
    x = [y for y in range(0, ref_length+1, bin_num)]
    for read in samfile.fetch(contig):
        rstart = read.reference_start
        rend = read.reference_end
        qstart = read.query_alignment_start
        qend = read.query_alignment_end
        qlength = read.infer_query_length()
        if qstart > buffer:
            trimmed_start = True
        if qend < qlength - buffer:
            trimmed_end = True
        if read.is_secondary:
            continue
        for i in range(max([0, rstart-bin_size//bin_step]), rend//bin_step+1):
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
            bases = min([rend, i*bin_step + bin_size]) - max([rstart, i*bin_step]) / bin_size
            if read.is_reverse:
                coverage_reverse += bases
            else:
                coverage_forward += bases
    return(coverage_forward, coverage_reverse, trimmed_starts, trimmed_ends, starts_in, ends_in, x)



def get_cov_stats_short(bamfile, contig):
    samfile = pysam.AlignmentFile(bamfile, 'rb')



def create_main_page(outfile, fasta, checkm_file, metabat_folder, long_bam, short_bam):
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
            bin_id, marker_lineage, genomes, markers, marker_sets, e0, e1, e2, e3, e4, e5, completeness, contamination, heterogeneity = line.split()
            checkm_dict[bin_id.split('.')[1]] = [marker_lineage, completeness, contamination, heterogeneity]
    outlist = []
    for i in os.listdir(metabat_folder):
        if not i.startswith('binned_contigs'):
            continue
        with open(os.path.join(metabat_folder, i)) as fa:
            bin = i.split('.')[1]
            ctgs = []
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
                if y >= bases_assembled/2
                    break
            n50 = j
            ctg_details = []
            bases_sequenced_long = 0
            bases_sequenced_short = 0
            for ctg in ctgs:
                cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in, ends_ind, x = get_cov_stats_long(long_bam, ctg)
                cov_forward, cov_reverse = get_cov_stats_short
                coverage_long = (sum(cov_forward) + sum(cov_reverse)) / len(cov_forward)
                bases_sequenced_long += coverage_long * len_dict[ctg]
                if short_bam is None:
                    coverage_ill = None
                    coverage_short = 0
                else:
                    coverage_ill = get_cov_stats_short(short_bam, ctg)
                    coverage_short = sum(coverage_ill) / len(coverage_ill)
                    bases_sequenced_long += coverage_short * len_dict[ctg]
                # create_contig_page(bin, ctg, cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in,
                #                    ends_ind, coverage_ill x)
                ctg_details.append((ctg,  len_dict[j], coverage_long, coverage_short))
            create_bin_page(bin, ctg_details)
            outlist.append([bin, max_contig, bases_assembled, n50, bases_sequenced_long/bases_assembled, bases_sequenced_short/bases_assembled] + checkm_dict[bin] + ['missing'])
    with open(outfile) as o:
        o.write(create_header(bin_list, 'index.html'))
        o.write(add_title("Overview of assembly", "assembled into " + len(bin_list) + " bins."))
        o.write(add_main("Overview of bins", "details for each bin", 12))
        o.write(create_table(['Bin', 'Max. contig (bp)', 'bases assembled', 'N50', 'average read depth (long)',
                              'average read depth (short)', 'marker lineage',  'completeness', 'Contamination',
                              'Heterozygosity', 'Best mash hit'], outlist))
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

    # generate short read page
    # generate long read page



def create_bin_page(page, fasta, bam):
    pass



checkm_file = snakemake.input.checkm_file
metabat_bins = snakemake.input.metabat_bins[:-5]
fasta = snakemake.input.fasta
try:
    os.makedirs('www')
except FileExistsError:
    pass


if os.path.exists("data/short_reads.fastq.gz"):
    short_bam = "data/final_cov.sort.bam"
    long_bam = "data/long_reads_for_web.bam"
else:
    long_bam = "data/final_cov.sort.bam"
    short_bam = None

create_main_page("www/index.html", fasta, checkm_file, metabat_folder, long_bam, short_bam)