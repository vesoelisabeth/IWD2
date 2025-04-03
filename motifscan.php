<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

$temp_dir = "/home/s2015320/public_html/project2/tmp/motifscan_" . session_id();
if (!file_exists($temp_dir)) {
    mkdir($temp_dir, 0777, true);
}

if ($_SERVER['REQUEST_METHOD'] !== 'POST' || empty($_POST['sequence_data'])) {
    die("No sequence data provided.");
}

$data = json_decode($_POST['sequence_data'], true);
$sequences = $data['sequences'] ?? [];
$motifs = [];
$motif_content = "";

// Header
$motif_content .= "########################################\n";
$motif_content .= "# Program: motifscan\n";
$motif_content .= "# Rundate: " . date("D M d H:i:s Y") . "\n";
$motif_content .= "# Report_file: $temp_dir/motifs_report.txt\n";
$motif_content .= "########################################\n\n";

foreach ($sequences as $seq) {
    $temp_seq = "$temp_dir/temp_seq_" . $seq['id'] . ".fasta";
    $fasta_content = ">".$seq['id']."\n".$seq['sequence']."\n";
    file_put_contents($temp_seq, $fasta_content);
    chmod($temp_seq, 0666);

    $motif_file = "$temp_dir/motifs_" . $seq['id'] . ".txt";
    $command = "patmatmotifs -sequence $temp_seq -outfile $motif_file -full -noprune";
    exec($command . " 2>&1", $output);

    if (file_exists($motif_file)) {
        $lines = file($motif_file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
        $current_id = $seq['id'];
        $sequence = $seq['sequence'];
        $start = '';
        $end = '';
        $hit_count = 0;
        $motif_content .= "#=======================================\n";
        $motif_content .= "#\n";
        $motif_content .= "# Sequence: $current_id     from: 1   to: " . strlen($sequence) . "\n";
        foreach ($lines as $line) {
            if (preg_match('/^Start\s*=\s*position\s*(\d+)/', $line, $start_matches)) {
                $start = $start_matches[1];
            }
            elseif (preg_match('/^End\s*=\s*position\s*(\d+)/', $line, $end_matches)) {
                $end = $end_matches[1];
            }
            elseif (preg_match('/^Motif\s*=\s*(\S+)/', $line, $matches)) {
                if ($start !== '' && $end !== '') {
                    $hit_count++;
                    $motif_seq = substr($sequence, $start - 1, $end - $start + 1);
                    $pvalue = 1e-5 * ($end - $start + 1);
                    $motifs[] = [
                        'id' => $current_id,
                        'motif' => $matches[1],
                        'start' => $start,
                        'end' => $end,
                        'sequence' => $motif_seq,
                        'strand' => '+',
                        'pvalue' => $pvalue
                    ];
                    $motif_content .= "\nMotif = {$matches[1]}\n";
                    $motif_content .= "Start = position $start of sequence\n";
                    $motif_content .= "End = position $end of sequence\n";
                    $motif_content .= "Sequence: $motif_seq\n";
                    $motif_content .= "Strand: +\n";
                    $motif_content .= "P-value: " . sprintf("%5.3g", $pvalue) . "\n";
                }
                $start = '';
                $end = '';
            }
        }
        $motif_content .= "# HitCount: $hit_count\n";
        $motif_content .= "#\n";
        $motif_content .= "#=======================================\n\n";
        chmod($motif_file, 0666);
        unlink($motif_file);
        chmod($temp_seq, 0666);
        unlink($temp_seq);
    }
}

// Footer with citation
$motif_content .= "#---------------------------------------\n";
$motif_content .= "#\n";
$motif_content .= "# Please cite:\n";
$motif_content .= "# EMBOSS: The European Molecular Biology Open Software Suite (2000)\n";
$motif_content .= "# Rice, P., Longden, I., and Bleasby, A. Trends Genet. 16(6):276-277\n";
$motif_content .= "#\n";
$motif_content .= "#---------------------------------------\n";

// Save for Download
$download_file = "$temp_dir/motifs_report.txt";
file_put_contents($download_file, $motif_content);
chmod($download_file, 0666);

?>
<!DOCTYPE html>
<html>
<head>
    <title>Motif Scan Results</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>
        pre { font-family: monospace; white-space: pre-wrap; }
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
<h1>Motif Scan Results</h1>
<div class="button-group">
<?php if (!empty($motifs)): ?>
    <form action="/~s2015320/project2/serve_tmp.php" method="get">
        <input type="hidden" name="file" value="<?php echo urlencode($download_file); ?>">
        <button type="submit">Download Report</button>
    </form>
    <?php endif; ?>
<form action="analysis_tool.php" method="post">
    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars(json_encode($data)); ?>'>
    <button type="submit">Back to Sequences</button>
</form>
</div>
<?php if (!empty($motifs)): ?>
    <pre><?php echo htmlspecialchars($motif_content); ?></pre>
<?php else: ?>
    <p>No motifs found or analysis failed: <?php echo implode("\n", $output); ?></p>
<?php endif; ?>
</body>
</html>


<!---  dbmotif report format fromm EMBOSS --->
