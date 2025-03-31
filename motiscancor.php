<?php 
session_start();

// Enable debugging
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// Debug log file
$debug_log = '/tmp/motif_debug.log';
file_put_contents($debug_log, "Script started\n", FILE_APPEND);

// Handle download request
$temp_output = "/home/s2015320/public_html/project2/motif_scan.out";
if (isset($_GET['download']) && $_GET['download'] == 'report') {
    if (file_exists($temp_output)) {
        header('Content-Type: text/plain');
        header('Content-Disposition: attachment; filename="motif_scan_report.txt"');
        header('Content-Length: ' . filesize($temp_output));
        readfile($temp_output);
        exit;
    } else {
        file_put_contents($debug_log, "Download failed: Report file not found\n", FILE_APPEND);
        die("Report file not found.");
    }
}

// Get sequence data from POST (only for scan)
if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    $sequence_data = json_decode($_POST['sequence_data'] ?? '{}', true);
    if (empty($sequence_data['sequences'])) {
        file_put_contents($debug_log, "No sequence data provided\n", FILE_APPEND);
        die("No sequence data provided.");
    }
    file_put_contents($debug_log, "Sequence data received\n", FILE_APPEND);

    // Store sequence data in session to preserve it
    $_SESSION['last_sequence_data'] = json_encode($sequence_data);

    // Prepare FASTA input for patmatmotifs
    $fasta_input = "";
    foreach ($sequence_data['sequences'] as $seq) {
        $fasta_input .= ">" . $seq['id'] . "\n" . $seq['sequence'] . "\n";
    }

    // Write FASTA to a fixed location
    $temp_fasta = "/home/s2015320/public_html/project2/temp_seq_" . uniqid() . ".fasta";
    file_put_contents($temp_fasta, $fasta_input);
    file_put_contents($debug_log, "FASTA Input:\n$fasta_input\n", FILE_APPEND);

    // Run patmatmotifs with -full and -noprune
    $prosite_path = "/usr/share/EMBOSS/data/PROSITE/prosite.lines";
    $command = "/usr/bin/patmatmotifs -sequence " . escapeshellarg($temp_fasta) . " -full -noprune -sformat fasta -rformat simple -outfile " . escapeshellarg($temp_output) . " 2>&1";
    $result = shell_exec($command);
    file_put_contents($debug_log, "Command: $command\nResult: $result\n", FILE_APPEND);

    if (!file_exists($temp_output) || filesize($temp_output) == 0) {
        $error = "Motif scan failed: " . htmlspecialchars($result);
        file_put_contents($debug_log, "Motif scan failed: $error\n", FILE_APPEND);
        unlink($temp_fasta);
        die($error);
    }
    file_put_contents($debug_log, "Motif scan output exists\n", FILE_APPEND);

    // Parse the simple output
    $output_content = file_get_contents($temp_output);
    $lines = explode("\n", $output_content);
    $motifs = [];
    $current_seq = '';
    $feature_data = [];

    foreach ($lines as $line) {
        if (preg_match('/^#\s+Sequence:\s+(\S+)/', $line, $matches)) {
            $current_seq = $matches[1];
        } elseif (preg_match('/^Feature:\s+(\d+)/', $line, $matches)) {
            $feature_data = []; // Reset for new feature
        } elseif (preg_match('/^Name:\s+(\S+)/', $line, $matches)) {
            $feature_data['name'] = $matches[1];
        } elseif (preg_match('/^Start:\s+(\d+)/', $line, $matches)) {
            $feature_data['start'] = $matches[1];
        } elseif (preg_match('/^End:\s+(\d+)/', $line, $matches)) {
            $feature_data['end'] = $matches[1];
        } elseif (preg_match('/^Motif:\s+(\S+)/', $line, $matches)) {
            $feature_data['motif'] = $matches[1];
            if (isset($feature_data['name'], $feature_data['start'], $feature_data['end'], $feature_data['motif'])) {
                $motifs[] = [
                    'sequence_id' => $current_seq,
                    'motif' => $feature_data['motif'],
                    'start' => $feature_data['start'],
                    'end' => $feature_data['end'],
                    'description' => $feature_data['name']
                ];
            }
        }
    }
    file_put_contents($debug_log, "Output: $output_content\nMotifs: " . print_r($motifs, true) . "\n", FILE_APPEND);

    // Clean up
   unlink($temp_fasta);
} else {

    // For GET requests (e.g., after download), load motifs from last run
    if (file_exists($temp_output)) {
        $output_content = file_get_contents($temp_output);
        $lines = explode("\n", $output_content);
        $motifs = [];
        $current_seq = '';
        $feature_data = [];

        foreach ($lines as $line) {
            if (preg_match('/^#\s+Sequence:\s+(\S+)/', $line, $matches)) {
                $current_seq = $matches[1];
            } elseif (preg_match('/^Feature:\s+(\d+)/', $line, $matches)) {
                $feature_data = []; // Reset for new feature
            } elseif (preg_match('/^Name:\s+(\S+)/', $line, $matches)) {
                $feature_data['name'] = $matches[1];
            } elseif (preg_match('/^Start:\s+(\d+)/', $line, $matches)) {
                $feature_data['start'] = $matches[1];
            } elseif (preg_match('/^End:\s+(\d+)/', $line, $matches)) {
                $feature_data['end'] = $matches[1];
            } elseif (preg_match('/^Motif:\s+(\S+)/', $line, $matches)) {
                $feature_data['motif'] = $matches[1];
                if (isset($feature_data['name'], $feature_data['start'], $feature_data['end'], $feature_data['motif'])) {
                    $motifs[] = [
                        'sequence_id' => $current_seq,
                        'motif' => $feature_data['motif'],
                        'start' => $feature_data['start'],
                        'end' => $feature_data['end'],
                        'description' => $feature_data['name']
                    ];
                }
            }
        }
        file_put_contents($debug_log, "Loaded motifs from file: " . print_r($motifs, true) . "\n", FILE_APPEND);
    } else {
        $motifs = [];
        file_put_contents($debug_log, "No previous motif output found\n", FILE_APPEND);
    }
}

file_put_contents($debug_log, "Reached HTML generation\n", FILE_APPEND);

// Generate HTML with navbar
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
    <div>
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
        <p>
            <a href="/~s2015320/project2/motifscan.php?download=report"><button>Download Report</button></a>
        </p>
        <form action="/~s2015320/project2/analysis_tool.php" method="post">
            <input type="hidden" name="sequence_data" value="<?php echo htmlspecialchars($_SESSION['last_sequence_data'] ?? ''); ?>">
            <button type="submit">Back to Sequences</button>
        </form>
    </div>
</body>
</html>

<!--  https://www.geeksforgeeks.org/how-to-trigger-a-file-download-when-clicking-an-html-button-or-javascript/  -->
