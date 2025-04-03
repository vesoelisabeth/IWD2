<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);
header('Content-Type: text/html; charset=UTF-8');

$tmp_dir = "/home/s2015320/public_html/project2/tmp";
$debug_log = "$tmp_dir/alignment_debug_" . session_id() . ".log";
if (!file_exists($tmp_dir)) mkdir($tmp_dir, 0777, true);

if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    file_put_contents($debug_log, "No sequence data\n", FILE_APPEND);
    die("No sequence data provided.");
}

$input_data = json_decode($_POST['sequence_data'], true);
$protein_name = $input_data['protein_name'] ?? 'Unknown Protein';
$taxonomy = $input_data['taxonomy'] ?? 'Unknown Taxonomy';

$temp_json = tempnam($tmp_dir, 'seq_');
file_put_contents($temp_json, $_POST['sequence_data']);
chmod($temp_json, 0666);
$script_path = "/home/s2015320/public_html/project2/clustal_analyze.py";
$python_path = "/home/s2015320/public_html/project2/myenv/bin/python3";
$command = "$python_path $script_path < $temp_json 2>&1";
exec($command, $output, $return_var);
$result = implode("\n", $output);
file_put_contents($debug_log, "Command: $command\nInput: {$_POST['sequence_data']}\nOutput: $result\nReturn: $return_var\n", FILE_APPEND);
unlink($temp_json);

if ($return_var !== 0 || empty($result)) {
    die("Alignment failed: " . htmlspecialchars($result));
}

$json_start = strpos($result, '{');
$result = $json_start !== false ? substr($result, $json_start) : '{}';
$data = json_decode($result, true);
if ($data === null) {
    die("JSON parsing failed: " . htmlspecialchars($result));
}

$conservation_score = $data['conservation_score'] ?? 'N/A';
$plot_url = basename($data['plot_url'] ?? ''); // Use basename
$msa_plot_url = basename($data['msa_plot_url'] ?? ''); // Use basename
$alignment = $data['alignment'] ?? [];

$aln_file = "$tmp_dir/alignment_" . session_id() . ".txt";
$aln_content = "Conservation Score: $conservation_score\n\n";
foreach ($alignment as $seq) {
    $aln_content .= "> {$seq['id']}\n{$seq['aligned_sequence']}\n";
}
file_put_contents($aln_file, $aln_content);
chmod($aln_file, 0666);

// --- Store Sequence Data for Back Button ---
$_SESSION['last_sequence_data'] = $_POST['sequence_data'];
?>
<!DOCTYPE html>
<html>
<head>
    <title>Conservation Analysis</title>
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
        <a href="credits.php"><button>Credits</button></a>
    </nav>
    <div>
        <h1>Conservation Analysis: <?php echo htmlspecialchars("$protein_name ($taxonomy)"); ?></h1>
        <div class="button-group">
            <form action="/~s2015320/project2/serve_tmp.php" method="get">
                <input type="hidden" name="file" value="<?php echo urlencode($aln_file); ?>">
                <button type="submit">Download Alignment</button>
            </form>
            <form method="POST" action="analysis_tool.php">
                <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>'>
                <button type="submit">Back</button>
            </form>
        </div>
        <h2>Conservation Score</h2>
        <p><?php echo htmlspecialchars($conservation_score); ?></p>
        <?php if ($plot_url && file_exists("$tmp_dir/$plot_url")): ?>
            <h2>Conservation Plot</h2>
            <img src="/~s2015320/project2/tmp/<?php echo htmlspecialchars($plot_url); ?>" alt="Conservation Plot" style="max-width: 100%;">
        <?php else: ?>
            <h2>Conservation Plot</h2>
            <p>Not generated: <?php echo htmlspecialchars(file_get_contents($debug_log)); ?></p>
        <?php endif; ?>
        <?php if ($msa_plot_url && file_exists("$tmp_dir/$msa_plot_url") && !empty($alignment)): ?>
            <h2>Multiple Sequence Alignment Plot</h2>
            <img src="/~s2015320/project2/tmp/<?php echo htmlspecialchars($msa_plot_url); ?>" alt="MSA Plot" style="max-width: 100%;">
        <?php else: ?>
            <h2>Multiple Sequence Alignment Plot</h2>
            <p>Not generated: <?php echo htmlspecialchars(file_get_contents($debug_log)); ?></p>
        <?php endif; ?>
        <?php if (!empty($alignment)): ?>
            <h2>Multiple Sequence Alignment</h2>
            <table border="1">
                <tr><th>ID</th><th>Aligned Sequence</th></tr>
                <?php foreach ($alignment as $seq): ?>
                    <tr>
                        <td><?php echo htmlspecialchars($seq['id']); ?></td>
                        <td><pre><?php echo htmlspecialchars($seq['aligned_sequence']); ?></pre></td>
                    </tr>
                <?php endforeach; ?>
            </table>
        <?php endif; ?>
    </div>
</body>
</html>
