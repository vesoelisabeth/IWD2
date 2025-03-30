<?php
session_start();

// Get sequence data from POST
$sequence_data = json_decode($_POST['sequence_data'] ?? '{}', true);
if (empty($sequence_data['sequences'])) {
    die("No sequence data provided.");
}

// Prepare FASTA input for patmatmotifs
$fasta_input = "";
foreach ($sequence_data['sequences'] as $seq) {
    $fasta_input .= ">" . $seq['id'] . "\n" . $seq['sequence'] . "\n";
}

// Write FASTA to a temporary file
$temp_fasta = tempnam(sys_get_temp_dir(), 'seq_') . '.fasta';
file_put_contents($temp_fasta, $fasta_input);

// Run patmatmotifs
$prosite_path = "/usr/local/share/EMBOSS/data/PROSITE/prosite.dat"; 
$temp_output = "/home/s2015320/public_html/project2/motif_scan.out";
$command = "patmatmotifs -sequence " . escapeshellarg($temp_fasta) ." -prune -sformat fasta -rformat simple -outfile " . escapeshellarg($temp_output) . " 2>&1";
$result = shell_exec($command);

// Debug: Log command and result
file_put_contents('/tmp/motif_debug.log', "Command: $command\nResult: $result\n", FILE_APPEND);

if (!file_exists($temp_output) || filesize($temp_output) == 0) {
    $error = "Motif scan failed: " . htmlspecialchars($result);
    unlink($temp_fasta);
    die($error);
}

// Parse the simple output
$output_content = file_get_contents($temp_output);
$lines = explode("\n", $output_content);
$motifs = [];
$current_seq = '';
foreach ($lines as $line) {
    if (preg_match('/^#\s+Sequence:\s+(\S+)/', $line, $matches)) {
        $current_seq = $matches[1];
    } elseif (preg_match('/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/', $line, $matches)) {
        $motifs[] = [
            'sequence_id' => $current_seq,
            'motif' => $matches[1],
            'start' => $matches[2],
            'end' => $matches[3],
            'description' => $matches[4]
        ];
    }
}

// Debug: Log parsed output
file_put_contents('/tmp/motif_debug.log', "Output: $output_content\nMotifs: " . print_r($motifs, true) . "\n", FILE_APPEND);

// Clean up
unlink($temp_fasta);
//unlink($temp_output);

// Generate HTML table with results
?>
<!DOCTYPE html>
<html>
<head>
    <title>Motif Scan Results</title>
    <link rel="stylesheet" href="styles_new.css">
<nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
     <!-- <a href="about.php"><button>About</button></a> --> 
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>

</head>
<body>
    <h1>Motif Scan Results</h1>
    <?php if (empty($motifs)): ?>
        <p>No motifs found in the provided sequences.</p>
    <?php else: ?>
        <table border="1">
            <tr>
                <th>Sequence ID</th>
                <th>Motif</th>
                <th>Start</th>
                <th>End</th>
                <th>Description</th>
            </tr>
            <?php foreach ($motifs as $motif): ?>
                <tr>
                    <td><?php echo htmlspecialchars($motif['sequence_id']); ?></td>
                    <td><?php echo htmlspecialchars($motif['motif']); ?></td>
                    <td><?php echo htmlspecialchars($motif['start']); ?></td>
                    <td><?php echo htmlspecialchars($motif['end']); ?></td>
                    <td><?php echo htmlspecialchars($motif['description']); ?></td>
                </tr>
            <?php endforeach; ?>
        </table>
    <?php endif; ?>
    <a href="/~s2015320/project2/analysis_tool.php"><button>Back to Query</button></a>
</body>
</html>
