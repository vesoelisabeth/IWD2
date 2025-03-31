<?php
// Enable debugging
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// Set proper headers
header('Content-Type: text/html; charset=UTF-8');

session_start();

// Debug log
$debug_log = '/home/s2015320/public_html/project2/alignment_debug.log';
file_put_contents($debug_log, "Script started\n", FILE_APPEND);

// Validate input
if ($_SERVER['REQUEST_METHOD'] !== 'POST') {
    file_put_contents($debug_log, "Invalid request method\n", FILE_APPEND);
    die("Only POST requests are allowed");
}

// Get sequence data from POST
$sequence_data = $_POST['sequence_data'] ?? '';
if (empty($sequence_data)) {
    file_put_contents($debug_log, "No sequence data provided\n", FILE_APPEND);
    die("No sequence data provided.");
}
file_put_contents($debug_log, "Sequence data: $sequence_data\n", FILE_APPEND);

// Write JSON to a temporary file
$temp_json = "/home/s2015320/public_html/project2/seq_" . uniqid() . ".json";
file_put_contents($temp_json, $sequence_data);
chmod($temp_json, 0644);
file_put_contents($debug_log, "Wrote JSON to $temp_json\n", FILE_APPEND);

// Set plot path
$plot_filename = "conservation_plot_" . session_id() . ".png";
$plot_path = "/home/s2015320/public_html/project2/" . $plot_filename;
$plot_url = "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/" . $plot_filename;
putenv("PLOT_PATH=$plot_path");


// Run Clustal Omega via Python script
$script_path = "/home/s2015320/public_html/project2/clustal_analyze.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "$python_path $script_path < $temp_json 2>&1";
exec($command, $output, $return_var);
$result = implode("\n", $output);
$result = preg_replace('/^Content-type: application\/json\s*/i', '', $result);
file_put_contents($debug_log, "Command: $command\nOutput: $result\nReturn: $return_var\n", FILE_APPEND);

// Clean up
unlink($temp_json);

// Check result
$data = json_decode($result, true);
if (!$data || $return_var !== 0) {
    file_put_contents($debug_log, "Alignment failed: $result\n", FILE_APPEND);
    die("Alignment failed: " . htmlspecialchars($result));
}

// Extract results
$conservation_score = $data['conservation_score'] ?? 'No score available';
$plot_url = $data['plot_url'] ?? '';
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
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
        <a href="about.php"><button>About</button></a>
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
            <p>No plot available (single sequence).</p>
        <?php endif; ?>
        <form action="analysis_tool.php" method="post">
            <input type="hidden" name="sequence_data" value="<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>">
            <button type="submit">Back to Sequences</button>
        </form>
    </div>
</body>
</html>
