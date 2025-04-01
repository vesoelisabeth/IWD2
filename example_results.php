<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

$job_id = "example_g6pase_aves";
$response = file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
    'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => json_encode(['action' => 'get_job', 'job_id' => $job_id])]
]));
$data = json_decode($response, true);
$files = $data['files'] ?? [];
?>
<!DOCTYPE html>
<html>
<head>
    <title>Example Results - Glucose-6-Phosphatase (Aves)</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>iframe { border: 1px solid #ccc; margin: 10px 0; }</style>
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
    <h2>Multiple Sequence Alignment (Conservation)</h2>
    <iframe src="https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/output/<?php echo basename($files['msa']); ?>" width="100%" height="400"></iframe>
    <h2>Motif Scan (PatMatMotif)</h2>
    <iframe src="https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/output/<?php echo basename($files['motifs']); ?>" width="100%" height="400"></iframe>
    <h2>Advanced Analysis</h2>
    <iframe src="https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/output/<?php echo basename($files['report']); ?>" width="100%" height="600"></iframe>
    <form action="analysis_tool.php" method="get">
        <button type="submit">Back to Analysis Tool</button>
    </form>
</body>
</html>
