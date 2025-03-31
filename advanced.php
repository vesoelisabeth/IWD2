<?php
session_start();

// --- Setup and Debugging ---
ini_set('display_errors', 1);  // Show errors on page
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// --- Input Validation ---
if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    die("No sequence data provided.");  // Exit if no data from analysis_tool.php
}

// --- Prepare Temporary Files ---
$temp_dir = "/tmp/advanced_" . session_id() . "/";  // Unique temp directory per session
mkdir($temp_dir, 0777, true);  // Create directory with full permissions
$temp_aln = $temp_dir . "input.aln";  // Alignment file path
$data = json_decode($_POST['sequence_data'], true);  // Decode incoming JSON data
$sequences = $data['sequences'] ?? [];  // Extract sequences
$fasta = "";
foreach ($sequences as $seq) {
    $fasta .= ">{$seq['id']}\n{$seq['sequence']}\n";  // Build FASTA format
}
file_put_contents($temp_dir . "input.fasta", $fasta);  // Write FASTA file

// --- Run Clustal Omega ---
$clustalo_cmd = "clustalo -i {$temp_dir}input.fasta -o $temp_aln --force --outfmt=clustal";
exec($clustalo_cmd . " 2>&1", $clustal_output);  // Align sequences, capture errors

// --- Run Advanced Analysis ---
$script_path = "/home/s2015320/public_html/project2/advanced_analysis.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "$python_path $script_path $temp_aln 2>&1";
exec($command, $output);  // Execute Python script, capture output
$report_path = array_pop($output);  // Last line is report file path
$report_url = str_replace("/tmp/", "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/tmp/", $report_path);  // Web-accessible URL
?>
<!DOCTYPE html>
<html>
<head>
    <title>Advanced Protein Analysis</title>
</head>
<body>
    <!-- Navigation -->
    <nav>
        <a href="indexx.html"><button>Home</button></a>
	<a href="about.php">button>About</button></a>        
	<a href="analysis_tool.php"><button>Analysis Tool</button></a>
	<a href="help.php"><button>Help</button></a>
	<a href="credits.php"><button>Statement of Credits></button></a>
    </nav>

    <!-- Page Header -->
    <h1>Advanced Protein Analysis</h1>

    <!-- Display Results -->
    <?php if (file_exists($report_path)): ?>
        <iframe src="<?php echo htmlspecialchars($report_url); ?>" width="100%" height="600"></iframe>
    <?php else: ?>
        <p>Analysis failed: <?php echo implode("\n", $output); ?></p>
    <?php endif; ?>

    <!-- Back Button -->
    <form action="analysis_tool.php" method="post">
        <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>'>
        <button type="submit">Back to Sequences</button>
    </form>
</body>
</html>
