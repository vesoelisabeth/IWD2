<?php
session_start();

//Specifically check the session state

$logged_in = isset($_SESSION['session_id']);
$results = '';
$data = null;

if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    if (isset($_POST['start_session'])) {
        $name = $_POST['name'] ?? '';
        $surname = $_POST['surname'] ?? '';
        if ($name && $surname) {
            $json_input = json_encode(['action' => 'start_session', 'name' => $name, 'surname' => $surname]);
            $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
            ]));
            if ($response === false) {
                $results = "Error: Failed to connect to sqlconnect.php.";
            } else {
                $data = json_decode($response, true);
                if (isset($data['session_id']) && $data['session_id']) {
                    $_SESSION['session_id'] = $data['session_id'];
                    $logged_in = true; // Update $logged_in here
                    $results = "Session started successfully!";
                } else {
                    $results = "Error: No valid session_id received. Response: " . htmlspecialchars($response);
                }
            }
        } else {
            $results = "Error: Name and surname are required.";
        }
     } elseif (isset($_POST['logout'])) {
        // Handle logout
        $session_id = $_SESSION['session_id'];
        $json_input = json_encode(['action' => 'end_session', 'session_id' => $session_id]);
        $context = stream_context_create([
            'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
        ]);
        @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, $context);
        session_unset();
        session_destroy();
        header("Location: analysis_tool.php");
        exit;


    } elseif (isset($_POST['fetch_proteins']) && $logged_in) { //Only can proceed if logged in
        $protein = $_POST['protein'] ?? '';
        $taxonomy = $_POST['taxonomy'] ?? '';
        $session_id = $_SESSION['session_id'] ?? '';
        if ($protein && $taxonomy && $session_id) { //Now the user is logged in so they may proceed
            $json_input = json_encode(['action' => 'search_protein', 'protein_name' => $protein, 'taxonomy' => $taxonomy, 'session_id' => $session_id]);
            $response = file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
            ]));
            $data = json_decode($response, true);
            if (!isset($data['error'])) {
                $results = "<h2>{$data['result']}</h2>";
                foreach ($data['sequences'] as $seq) {
                    $results .= "<div class='sequence'><strong>{$seq['id']}</strong><br>{$seq['description']}<br><p>Organism: {$seq['organism']}</p><pre>{$seq['sequence']}</pre></div>";
                }
            } else {
                $results = "Error: {$data['error']}";
	   }
        } else {
            $results = "Error: Please provide protein and taxonomy.";
        }
    } elseif (isset($_POST['fetch_proteins']) && !$logged_in) {
        $results = "Error: Please log in first.";
    }
}
?>
<!DOCTYPE html>
<html>
<head>
    <title>Protein Fetcher - Analysis Tool</title>
    <link rel="stylesheet" href="styles_new.css">
<style>
        .header-container {
            display: flex;
            align-items: center;
            justify-content: flex-start;
            gap: 10px; /* Space between elements */
            margin-bottom: 20px; /* Space below the header */
        }
        .header-container button {
            padding: 8px 16px;
            font-size: 14px;
        }
    </style>
</head>
<body>
    <nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
     <!-- <a href="about.php"><button>About</button></a> --> 
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>
    <h1>Analysis Tool</h1>
    <?php if (!$logged_in): ?>
        <div id="login">
            <h2>Login</h2>
		<p> You need to provide your details for optimised user experience. </p>

            <form method="POST">
                <div class="input-group">
                    Name: <input name="name" type="text" placeholder="Enter your name"><br>
                </div>
                <div class="input-group">
                    Surname: <input name="surname" type="text" placeholder="Enter your surname"><br>
                </div>
                <button type="submit" name="start_session">Start</button>
            </form>
        </div>
    <?php else: ?>
      <div id="logout">
            <form method="POST">
                <button type="submit" name="logout">Log Out</button>
            </form>
        </div>
       <div id="query">
            <h2>Query Proteins</h2>
		<p> Insert your protein and taxonomy of interest and we will tell you more about it. </p>
            <form method="POST">
                <div class="input-group">
                     Protein: <input name="protein" type="text" placeholder="e.g., glucose-6-phosphatase"><br>
                </div>
                <div class="input-group">
                    Taxonomy: <input name="taxonomy" type="text" placeholder="e.g., Aves"><br>
                </div>
                <button type="submit" name="fetch_proteins">Fetch</button>
            </form>
        </div>
        <div id="results">
            <?php echo $results; ?>
            <?php if ($data && !empty($data['sequences'])): ?>
                <form method="POST" action="analyze.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars(json_encode($data)); ?>'>
                    <button type="submit">Analyze Conservation</button>
                </form>
                <form method="POST" action="motifscan.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars(json_encode($data)); ?>'>
                    <button type="submit">Scan Motifs</button>
                </form>
            <?php endif; ?>
        </div>
    <?php endif; ?>
</body>
</html>

<!--- for recording sessions: https://stackoverflow.com/questions/2823913/how-do-i-create-a-user-history -->
 
