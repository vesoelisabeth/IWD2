<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

$temp_dir = "/home/s2015320/public_html/project2/tmp/motifscan_" . session_id();
$debug_log = "$temp_dir/motifscan_debug.log";
mkdir($temp_dir, 0777, true);

if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    file_put_contents($debug_log, "No sequence data\n", FILE_APPEND);
    die("No sequence data provided.");
}

$data = json_decode($_POST['sequence_data'], true);
$sequences = $data['sequences'] ?? [];
$temp_seq = "$temp_dir/temp_seq.fasta";
$fasta_content = "";
foreach ($sequences as $seq) {
    $fasta_content .= ">".$seq['id']."\n".$seq['sequence']."\n";
}
file_put_contents($temp_seq, $fasta_content);
chmod($temp_seq, 0666);

// Run Motif Scan
$motif_file = "$temp_dir/motifs.txt";
$command = "patmatmotifs -sequence $temp_seq -outfile $motif_file -full -noprune";
exec($command . " 2>&1", $output);
$motif_output = file_exists($motif_file) ? file_get_contents($motif_file) : "Not found";
file_put_contents($debug_log, "Command: $command\nFasta: $fasta_content\nOutput: " . implode("\n", $output) . "\nMotif File:\n$motif_output\n", FILE_APPEND);

// Parse Results with Position
$motifs = [];
$motif_content = "Motif Scan Results\n\n";
if (file_exists($motif_file)) {
    $lines = file($motif_file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
    $current_id = '';
    foreach ($lines as $line) {
        if (preg_match('/^Sequence:\s*(\S+)/', $line, $id_matches)) {
            $current_id = $id_matches[1];
            continue;
        }
        if (preg_match('/^Motif\s*=\s*(\S+)\s*(?:\[(\d+)-(\d+)\])?/', $line, $matches)) {
            $motifs[] = [
                'id' => $current_id,
                'motif' => $matches[1],
                'start' => isset($matches[2]) ? $matches[2] : '',
                'end' => isset($matches[3]) ? $matches[3] : ''
            ];
            $motif_content .= "ID: $current_id, Motif: {$matches[1]}" . (isset($matches[2]) ? ", Start: {$matches[2]}, End: {$matches[3]}" : "") . "\n";
        }
    }
}

// Save for Download
$download_file = "$temp_dir/motifs_download.txt";
file_put_contents($download_file, $motif_content);
chmod($download_file, 0666);

?>
<!DOCTYPE html>
<html>
<head>
    <title>Motif Scan Results</title>
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
<h1>Motif Scan Results</h1>
<?php if (!empty($motifs)): ?>
    <table border="1">
        <tr><th>Sequence ID</th><th>Motif</th><th>Start</th><th>End</th></tr>
        <?php foreach ($motifs as $motif): ?>
            <tr>
                <td><?php echo htmlspecialchars($motif['id']); ?></td>
                <td><?php echo htmlspecialchars($motif['motif']); ?></td>
                <td><?php echo htmlspecialchars($motif['start']); ?></td>
                <td><?php echo htmlspecialchars($motif['end']); ?></td>
            </tr>
        <?php endforeach; ?>
    </table>
    <form action="/~s2015320/project2/serve_tmp.php" method="get">
        <input type="hidden" name="file" value="<?php echo urlencode(basename($download_file)); ?>">
        <button type="submit">Download Motifs</button>
    </form>
<?php else: ?>
    <p>No motifs found or analysis failed: <?php echo implode("\n", $output); ?></p>
<?php endif; ?>
<form action="analysis_tool.php" method="post">
    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>'>
    <button type="submit">Back to Sequences</button>
</form>
<?php
// Cleanup
$files = glob("$temp_dir/*");
foreach ($files as $file) {
    if (is_file($file) && $file !== $download_file) unlink($file);
}
rmdir($temp_dir);
?>
</body>
</html>
