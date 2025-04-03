<?php
session_start();

// --- Setup and Debugging ---
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);
ini_set('max_execution_time', 600); // 10 minutes timeout

// --- Input Validation ---
if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    die("No sequence data provided.");
}

// --- Prepare Output Directory ---
$temp_dir = "/home/s2015320/public_html/project2/tmp/advanced_" . session_id();
if (!file_exists($temp_dir)) {
    mkdir($temp_dir, 0777, true); // Session-specific, writable
}
$session_id = session_id();
$temp_fasta = "$temp_dir/input.fasta";
$temp_aln = "$temp_dir/input.aln";
$temp_dnd = "$temp_dir/temp.dnd"; // For phylogeny tree
$report_path = "$temp_dir/report.html";
$report_url = "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/tmp/advanced_$session_id/report.html";

$data = json_decode($_POST['sequence_data'], true);
$sequences = $data['sequences'] ?? [];
$fasta = "";
foreach ($sequences as $seq) {
    $fasta .= ">{$seq['id']}\n{$seq['sequence']}\n";
}
if (!file_put_contents($temp_fasta, $fasta)) {
    die("Failed to write FASTA file: $temp_fasta - Check permissions.");
}
chmod($temp_fasta, 0666);

// --- Run Clustal Omega with Tree Output ---
$clustalo_cmd = "clustalo -i $temp_fasta -o $temp_aln --force --outfmt=clustal --guidetree-out=$temp_dnd";
exec($clustalo_cmd . " 2>&1", $clustal_output);

// --- Run Advanced Analysis ---
$script_path = "/home/s2015320/public_html/project2/advanced_analysis.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "MPLCONFIGDIR=/tmp/matplotlib_cache SESSION_ID=$session_id TEMP_DIR=$temp_dir $python_path $script_path $temp_aln 2>&1";
exec($command, $output);

// --- Store Sequence Data for Back Button ---
$_SESSION['last_sequence_data'] = $_POST['sequence_data'];
?>
<!DOCTYPE html>
<html>
<head>
    <title>Advanced Protein Analysis</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>
        pre { font-family: monospace; white-space: pre-wrap; }
        .button-group { margin-bottom: 10px; }
        .button-group form { display: inline; margin-right: 10px; }
    </style>
</head>
<body>
<nav class="top-nav">
    <a href="indexx.html"><button>Home</button></a>
    <a href="about.php"><button>About</button></a>
    <a href="analysis_tool.php"><button>Analysis Tool</button></a>
    <a href="help.php"><button>Help</button></a>
    <a href="credits.php"><button>Statement of Credits</button></a>
</nav>
<h1>Advanced Protein Analysis</h1>
<div class="button-group">
    <?php if (file_exists($report_path)): ?>
        <form action="/~s2015320/project2/serve_tmp.php" method="get">
            <input type="hidden" name="file" value="<?php echo urlencode($report_path); ?>">
            <button type="submit">Download Report</button>
        </form>
    <?php endif; ?>
    <form action="analysis_tool.php" method="post">
        <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>'>
        <button type="submit">Back to Sequences</button>
    </form>
</div>
<?php if (file_exists($report_path)): ?>
    <iframe src="<?php echo htmlspecialchars($report_url); ?>" width="100%" height="600"></iframe>
<?php else: ?>
    <p>Analysis failed: <pre><?php echo htmlspecialchars(implode("\n", $output)); ?></pre></p>
<?php endif; ?>
</body>
</html>
