<?php
session_start();

// Enable debugging
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// Set proper headers
header('Content-Type: text/html; charset=UTF-8');


// Debug log (defined before use)
$debug_log = '/tmp/alignment_debug_' . session_id() . '.log';


// Validate input
if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    if (file_exists('/tmp')) {
        file_put_contents($debug_log, "Invalid or no sequence data\n", FILE_APPEND);
    }
    die("No sequence data provided.");
}

// Temp JSON in /tmp
$temp_json = tempnam('/tmp', 'seq_');
file_put_contents($temp_json, $_POST['sequence_data']);

// Set plot paths
$plot_filename = "conservation_plot_" . session_id() . ".png";
$plot_path = "/home/s2015320/public_html/project2/" . $plot_filename;
$plot_url = "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/" . $plot_filename;
if (!file_exists($plot_path)) {
    @touch($plot_path); // Try to create, ignore failure
}
if (!is_writable($plot_path)) {
    $plot_path = "/tmp/" . $plot_filename; // Fallback to /tmp
}
putenv("PLOT_PATH=$plot_path");

$msa_plot_filename = "msa_plot_" . session_id() . ".png";
$msa_plot_path = "/home/s2015320/public_html/project2/" . $msa_plot_filename;
if (!file_exists($msa_plot_path)) {
    @touch($msa_plot_path); // Try to create
}
if (!is_writable($msa_plot_path)) {
    $msa_plot_path = "/tmp/" . $msa_plot_filename; // Fallback
}

// Run Clustal Omega via Python script
$script_path = "/home/s2015320/public_html/project2/clustal_analyze.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "$python_path $script_path < $temp_json 2>/dev/null";
exec($command, $output, $return_var);
$result = implode("\n", $output);
$result = preg_replace('/^Content-type: application\/json\s*/i', '', $result);


// Clean up
unlink($temp_json);

// Check result
$data = json_decode($result, true);
if (!$data || $return_var !== 0) {
    file_put_contents($debug_log, "Alignment failed: $result\n", FILE_APPEND);
  }
    die("Alignment failed: " . htmlspecialchars($result));

// Extract results
$conservation_score = $data['conservation_score'] ?? 'No score available';
$plot_url = $data['plot_url'] ?? $plot_url;
$msa_plot_url = $data['msa_plot_url'] ?? "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/" . $msa_plot_filename;
$alignment = $data['alignment'] ?? [];

// Output results
?>
<!DOCTYPE html>
<html>
<head>
    <title>Conservation Analysis</title>
    <link rel="stylesheet" href="styles_new.css">
</head>
<body>
<nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
	<a href="about.php"><button>About</button></a>        
	<a href="analysis_tool.php"><button>Analysis Tool</button></a>
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>
    <div>
<h1>Conservation Analysis</h1>
        <h2>Conservation Score</h2>
        <p><?php echo htmlspecialchars($conservation_score); ?></p>
        <?php if ($plot_url): ?>
            <h2>Conservation Plot</h2>
            <img src="<?php echo htmlspecialchars($plot_url); ?>" alt="Conservation Plot">
        <?php else: ?>
            <p>No conservation plot available (single sequence).</p>
        <?php endif; ?>
        <?php if ($msa_plot_url && !empty($alignment)): ?>
            <h2>Multiple Sequence Alignment Plot</h2>
            <img src="<?php echo htmlspecialchars($msa_plot_url); ?>" alt="MSA Plot">
        <?php endif; ?>
        <?php if (!empty($alignment)): ?>
            <h2>Multiple Sequence Alignment</h2>
            <table border="1">
                <tr><th>Sequence ID</th><th>Aligned Sequence</th></tr>
                <?php foreach ($alignment as $seq): ?>
                    <tr>
                        <td><?php echo htmlspecialchars($seq['id']); ?></td>
                        <td><pre><?php echo htmlspecialchars($seq['aligned_sequence']); ?></pre></td>
                    </tr>
<?php endforeach; ?>
            </table>
        <?php endif; ?>
        <form action="analysis_tool.php" method="post">
            <input type="hidden" name="sequence_data" value="<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>">
            <button type="submit">Back to Sequences</button>
        </form>
    </div>
</body>
</html>
