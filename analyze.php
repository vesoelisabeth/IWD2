<?php
// Enable debugging
ini_set('display_errors', 1);
error_reporting(E_ALL);

// Set proper headers
header('Content-Type: text/html; charset=UTF-8');

session_start();

// Validate input
if ($_SERVER['REQUEST_METHOD'] !== 'POST') {
    die("Only POST requests are allowed");
}

// Get sequence data from POST
$sequence_data = $_POST['sequence_data'] ?? '';
if (empty($sequence_data)) {
    die("No sequence data provided.");
}

// Write JSON to a temporary file
$temp_json = tempnam(sys_get_temp_dir(), 'seq_') . '.json';
file_put_contents($temp_json, $sequence_data);

// Run Clustal Omega via Python script
$script_path = "/home/s2015320/public_html/project2/clustal_analyze.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "$python_path $script_path < $temp_json 2>/dev/null";
exec($command, $output, $return_var);
$result = implode("\n", $output);
$result = preg_replace('/^Content-type: application\/json\s*/i', '', $result);
$data = json_decode($result, true);

// Clean up
unlink($temp_json);

// Check result
if (!$data) {
    die("Alignment failed: " . htmlspecialchars($result));
}

// Extract results
$conservation_score = $data['conservation_score'] ?? 'No score available';
$plot_url = $data['plot_url'] ?? '';

// Output results
?>
<!DOCTYPE html>
<html>
<head>
    <title>Conservation Analysis</title>
    <link rel="stylesheet" href="styles_new.css">
</head>
<body>
    <h1>Conservation Analysis</h1>
    <h2>Conservation Score</h2>
    <p><?php echo htmlspecialchars($conservation_score); ?></p>
    <?php if ($plot_url): ?>
        <h2>Conservation Plot</h2>
        <img src="<?php echo htmlspecialchars($plot_url); ?>" alt="Conservation Plot">
    <?php else: ?>
        <p>No plot available (single sequence).</p>
    <?php endif; ?>
    <a href="analysis_tool.php"><button>Back to Query</button></a>
</body>
</html>
