<?php
session_start();

// --- Setup and Debugging ---
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    die("No sequence data provided.");
}

// --- Prepare Temporary Files ---
$temp_dir = "/tmp/motifscan_" . session_id() . "/";
mkdir($temp_dir, 0777, true);
$temp_seq = $temp_dir . "temp_seq.fasta";
$data = json_decode($_POST['sequence_data'], true);
$sequences = $data['sequences'] ?? [];
$fasta_content = "";
foreach ($sequences as $seq) {
    $fasta_content .= ">{$seq['id']}\n{$seq['sequence']}\n";
}
file_put_contents($temp_seq, $fasta_content);

// --- Run Motif Scan ---
$command = "patmatmotifs -sequence $temp_seq -outfile $temp_dir/motifs.txt";
exec($command . " 2>&1", $output);

// --- Parse Results ---
$motif_file = $temp_dir . "motifs.txt";
$motifs = [];
if (file_exists($motif_file)) {
    $lines = file($motif_file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
    foreach ($lines as $line) {
        if (preg_match('/^Motif\s*=\s*(\S+)/', $line, $matches)) {
            $motifs[] = $matches[1];
        }
    }
}

// --- Generate HTML ---
?>
<!DOCTYPE html>
<html>
<head>
    <title>Motif Scan Results</title>
    <link rel="stylesheet" href="styles_new.css">
</head>
<body>
    <nav>
        <a href="indexx.html"><button>Home</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
    </nav>
    <h1>Motif Scan Results</h1>
    <?php if (!empty($motifs)): ?>
        <table border="1">
            <tr><th>Motif</th></tr>
            <?php foreach ($motifs as $motif): ?>
                <tr><td><?php echo htmlspecialchars($motif); ?></td></tr>
            <?php endforeach; ?>
        </table>
    <?php else: ?>
        <p>No motifs found or analysis failed: <?php echo implode("\n", $output); ?></p>
    <?php endif; ?>
    <form action="analysis_tool.php" method="post">
        <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>'>
        <button type="submit">Back to Sequences</button>
    </form>
    <?php
    // --- Cleanup ---
    unlink($temp_seq);
    unlink($motif_file);
    rmdir($temp_dir);
    ?>
</body>
</html>
