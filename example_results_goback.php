<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

$job_id = "example_g6pase_aves";
$base_url = "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/output/";
$json_path = "/home/s2015320/public_html/project2/output/conservation_{$job_id}.json";
$json_content = file_get_contents($json_path);


$json_content = preg_replace('/^Content-type: application\/json\s*/', '', $json_content);
$json_data = json_decode($json_content, true);
?>
<!DOCTYPE html>
<html>
<head>
    <title>Example Results - Glucose-6-Phosphatase (Aves)</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>iframe, img { border: 1px solid #ccc; margin: 10px 0; max-width: 100%; }</style>
</head>
<body>
    <nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
        <a href="about.php"><button>About</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>
    <h1>Example Results: Glucose-6-Phosphatase (Aves)</h1>
    <h2>Conservation Analysis</h2>
    <p>Conservation Score: <?php echo htmlspecialchars($json_data['conservation_score'] ?? 'N/A'); ?></p>
    <img src="<?php echo $base_url . "conservation_plot_{$job_id}.png"; ?>" alt="Conservation Plot">
    <img src="<?php echo $base_url . "msa_plot_{$job_id}.png"; ?>" alt="MSA Plot">
    <h2>Multiple Sequence Alignment</h2>
    <iframe src="<?php echo $base_url . "msa_{$job_id}.aln"; ?>" width="100%" height="400"></iframe>
    <h2>Motif Scan (PatMatMotif)</h2>
    <iframe src="<?php echo $base_url . "motifs_{$job_id}.txt"; ?>" width="100%" height="400"></iframe>
    <h2>Advanced Analysis</h2>
    <iframe src="<?php echo $base_url . "report.html"; ?>" width="100%" height="600"></iframe>
    <form action="analysis_tool.php" method="get">
        <button type="submit">Back to Analysis Tool</button>
    </form>
</body>
</html>
