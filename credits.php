<!DOCTYPE html>
<html>
<head>
    <title>Credits - Protein Fetcher Tool</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>
        body { 
            font-family: Arial, sans-serif; 
            line-height: 1.6; 
            width: 90%; 
            max-width: 1200px; 
            margin: 0 auto; 
            padding: 20px; 
            min-height: 100vh; 
            display: flex; 
            flex-direction: column; 
            justify-content: space-between; 
        }
        h1 { 
            color: #006400; 
            text-align: center; 
            margin-bottom: 30px; 
        }
        h2 { 
            color: #006400; 
            margin-top: 40px; 
            border-bottom: 2px solid #ecf0f1; 
            padding-bottom: 8px; 
        }
        p { 
            margin: 15px 0; 
        }
        ul { 
            margin: 15px 0 25px 0; 
            padding-left: 30px; 
        }
        li { 
            margin-bottom: 15px; 
        }
        a { 
            color: #006400; 
            text-decoration: underline; 
        }
        .top-nav { 
            text-align: center; 
            margin-bottom: 30px; 
            width: 100vw; 
            position: relative; 
            left: 50%; 
            transform: translateX(-50%); 
            padding: 10px 0; 
        }
        footer { 
            text-align: center; 
            padding: 10px 0; 
            border-top: 1px solid #ecf0f1; 
            width: 100%; 
        }
    </style>
</head>
<body>
    <nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
        <a href="about.php"><button>About</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>
    <h1>Statement of Credits</h1>

    <h2>Sources of Code</h2>
    <p>The Protein Fetcher Tool was built on various open resources, adapted for its specific needs:</p>
    <ul>
        <li><a href="https://www.w3schools.com/php/php_sessions.asp">W3Schools PHP Sessions</a> (Accessed on: 03 April 2025) provided foundational session management techniques, used in `analysis_tool.php` and `analyze.php` to track user inputs across pages.</li>
        <li><a href="https://www.geeksforgeeks.org/how-to-call-python-script-from-php/">GeeksforGeeks PHP-Python Integration</a> (Accessed on: 03 April 2025) guided subprocess calls from PHP to Python, implemented in `analyze.php`, `motifscan.php`, and `advanced.php` for running `clustal_analyze.py` and `advanced_analysis.py`.</li>
        <li><a href="https://stackoverflow.com/questions/39306676/how-to-execute-a-python-script-from-php">Stack Overflow Python Execution</a> (Accessed on: 03 April 2025) offered solutions for executing Python scripts with environment variables, applied in `generate_example.py` and `biopython_connect.py` for sequence retrieval.</li>
        <li><a href="https://www.w3schools.com/php/php_file_upload.asp">W3Schools File Handling</a> (Accessed on: 03 April 2025) informed temporary file management, used in `sqlconnect.php` and `generate_example.py` for storing and retrieving analysis outputs.</li>
        <li><a href="https://www.geeksforgeeks.org/how-to-create-a-dynamic-json-file-by-fetching-data-from-localstorage-using-php/">GeeksforGeeks JSON Handling</a> (Accessed on: 03 April 2025) shaped JSON output generation, adapted in `clustal_analyze.py` for conservation results.</li>
        <li><a href="https://www.w3schools.com/css/css_navbar_horizontal.asp">W3Schools CSS Horizontal Navbar</a> (Accessed on: 03 April 2025) provided the structure for the horizontal navigation bar, implemented in `indexx.html` for site-wide navigation.</li>
        <li><a href="https://www.w3schools.com/html/html_favicon.asp">W3Schools HTML Favicon</a> (Accessed on: 03 April 2025) guided favicon integration, used in `indexx.html` to enhance the site’s branding.</li>
        <li><a href="https://biopython.org/docs/latest/Tutorial/chapter_msa.html">BioPython MSA Tutorial</a> (Accessed on: 03 April 2025) informed multiple sequence alignment techniques, adapted in `advanced_analysis.py` for MSA and phylogenetic tree generation.</li>
        <li><a href="https://www.biostars.org/p/382859/">BioStars Clustal Omega Post</a> (Accessed on: 03 April 2025) served as a code source for files using Clustal Omega analysis, such as `clustal_analyze.py`, aiding in structuring multiple sequence alignment workflows.</li>
        <li><a href="https://emboss.bioinformatics.nl/cgi-bin/emboss/help/patmatmotifs">EMBOSS PatMatMotifs Help</a> (Accessed on: 03 April 2025) provided guidance for motif scanning, implemented in `motifscan.php` and `generate_example.py` for identifying protein motifs.</li>
        <li><a href="https://emboss.sourceforge.net/docs/themes/ReportFormats.html">EMBOSS Report Formats</a> (Accessed on: 03 April 2025) informed report formatting for motif analysis outputs, used in `motifscan.php`.</li>
        <li><a href="https://biopython.org/docs/latest/Tutorial/chapter_motifs.html">BioPython Motifs Tutorial</a> (Accessed on: 03 April 2025) provided techniques for motif identification, adapted in `motifscan.php` for motif scanning workflows.</li>
        <li><a href="https://www.w3schools.com/python/matplotlib_pyplot.asp">W3Schools Matplotlib Pyplot</a> (Accessed on: 03 April 2025) and associated pages offered plotting techniques, used in `clustal_analyze.py`, `advanced_analysis.py`, and `generate_example.py` for generating conservation and phylogenetic plots.</li>
        <li><a href="https://www.w3schools.com/php/php_mysql_connect.asp">W3Schools PHP MySQL Connect</a> (Accessed on: 03 April 2025) provided MySQL connection basics, used in `sqlconnect.php` for database integration, alongside BPSM Biopython Lecture Notes for NCBI retrieval concepts, BPSM/IWD2 SQL Lectures for query structuring, and <a href="https://stackoverflow.com/questions/52535136/python-return-json-to-php-and-use-json-decode">Stack Overflow JSON to PHP</a> (Accessed on: 03 April 2025) for JSON decoding techniques, all influencing `sqlconnect.php`.</li>
        <li><a href="https://github.com/biopython/biopython">BioPython GitHub Repository</a> (Accessed on: 03 April 2025) provided the core BioPython libraries and examples, influencing `biopython_connect.py`, `clustal_analyze.py`, and `advanced_analysis.py` for sequence retrieval and analysis.</li>
        <li><a href="https://github.com/Ensembl/ensembl-webcode">Ensembl Webcode GitHub</a> (Accessed on: 03 April 2025) offered inspiration for PHP-MySQL-BioPython integration, influencing `analysis_tool.php` and `sqlconnect.php` for web-based bioinformatics workflows.</li>
        <li><a href="https://github.com/galaxyproject/galaxy">Galaxy Project GitHub</a> (Accessed on: 03 April 2025) provided examples of Python-based bioinformatics tools with web interfaces, adapted in `analysis_tool.php` and `advanced_analysis.py` for tool orchestration and visualization.</li>
        <li>BPSM and IWD2 Lecture Notes (University of Edinburgh) supplied foundational concepts for sequence analysis and database integration, influencing `biopython_connect.py` for NCBI retrieval and `sqlconnect.php` for MySQL queries (BPSM), while providing basic code and theory for front- and backend structuring in files like `analysis_tool.php` and `advanced_analysis.py` (IWD2).</li>
    </ul>

    <h2>AI Tools Used</h2>
    <p>The development leveraged AI assistance for specific tasks:</p>
    <ul>
        <li><strong>Grok (xAI)</strong> was used for extensive debugging, such as resolving syntax errors in `analyze.php` (e.g., unexpected token issues), fixing subprocess failures in `generate_example.py`, and correcting file path issues in `clustal_analyze.py`. It created error traps, like retry logic for NCBI API calls in `biopython_connect.py`, adjusted basic form handling in `analysis_tool.php` for user inputs, and expanded code for complex tasks, such as adding Shannon entropy calculations in `clustal_analyze.py` and phylogenetic tree generation in `advanced_analysis.py`. Additionally, Grok suggested using the `tmp` directory format in `generate_example.py` and `sqlconnect.php` to bypass server permission restrictions for file storage, assisted with syncing JSON outputs, AJAX calls, and formatting the history table in `analysis_tool.php` to ensure seamless data flow and user interaction, helped implement `basename` in `clustal_analyze.py` and `analysis_tool.php` when troubleshooting the lack of fluid transfer between plot display from backend to frontend, and generated and modified `styles_new.css` to refine the site’s visual styling and layout consistency.</li>
    </ul>

    <h2>Scientific References</h2>
    <p>The tool’s biological insights draw from:</p>
    <ul>
        <li>Sievers, F. and Higgins, D.G., 2018. Clustal Omega for making accurate alignments of many protein sequences. *Protein Science*, 27(1), pp.135-145—basis for MSA in conservation analysis.</li>
        <li>Bywater, R.P., 2015. Prediction of protein structural features from sequence data based on Shannon entropy and Kolmogorov complexity. *PloS One*, 10(4), p.e0119306—foundation for entropy-based conservation scoring.</li>
        <li>Mithieux, G., Rajas, F. and Gautier-Stein, A., 2014. Glucose-6-phosphatase and glucose homeostasis: From metabolic control to therapeutic perspectives. *Diabetes*, 63(12), pp.3995-4006—context for the example dataset.</li>
    </ul>

    <footer>
        <p>For help, contact [your email].</p>
    </footer>
</body>
</html>
