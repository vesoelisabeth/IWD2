<!DOCTYPE html>
<html>
<head>
    <title>Help - Protein Fetcher Tool</title>
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
        ul { 
            margin: 15px 0 25px 0; 
            padding-left: 30px; 
        }
        li { 
            margin-bottom: 15px; 
        }
        p { 
            margin: 15px 0; 
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
        .example-image { 
            width: 200px; 
            height: auto; 
            display: block; 
            margin: 10px auto; 
        }
        footer { 
            text-align: center; 
            margin-top: 80px; 
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
    <h1>Help & Context: Protein Fetcher Tool</h1>

    <h2>Overview</h2>
    <p>The Protein Fetcher Tool assists biologists by analyzing protein sequences from the NCBI database, providing insights into conservation, motifs, phylogeny, and structure to support research discoveries.</p>

    <h2>Sequence Retrieval</h2>
    <p>The tool employs BLAST to fetch sequences by name and taxonomy, revealing protein roles across species while supporting studies of evolutionary diversity and functional variations in metabolism or disease.</p>

    <h2>Conservation Analysis</h2>
    <p>Using Clustal Omega for multiple sequence alignment (MSA), the tool identifies conserved regions critical for function (Sievers & Higgins, 2018). Additionally, Shannon entropy scores and plots map variability, helping pinpoint key residues for further study (Bywater, 2015).</p>

    <h2>Motif Identification</h2>
    <p>It detects motifs like phosphorylation sites, which indicate regulatory or structural roles, and highlights these functional elements to aid in pathway analysis.</p>

    <h2>Advanced Insights</h2>
    <p>The tool generates phylogenetic trees from MSA data to trace protein evolution across species, predicts folding patterns to study stability and interactions, and provides PDB models to visualise active sites and binding pockets.</p>

    <h2>Example Dataset</h2>
    <p>The pre-analysed glucose-6-phosphatase dataset from Aves offers a preview of the tool’s outputs. It includes glucose-6-phosphate, its substrate—explore it <a href="https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/analysis_tool.php">here</a> or learn more on <a href="https://en.wikipedia.org/wiki/Glucose_6-phosphatase">Wikipedia</a>. <img class="example-image" src="Alpha-D-glucopyranose_6-phosphate.svg (1).png" alt="Glucose-6-phosphate"></p>

    <h2>Applications</h2>
    <p>It explores evolution, function, and structure for regulation or disease research, guiding experiments targeting conserved regions or orthologs.</p>

    <footer>
        <p>For help, contact [s2015320@ed.ac.uk]. Log in, enter protein and taxonomy, select options, and Enjoy the tool!</p>
    </footer>
</body>
</html>
