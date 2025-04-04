<?php
session_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

// Step 1: Initialize variables
$logged_in = isset($_SESSION['session_id']);
$results = '';
$history = [];
$debug_log = "/tmp/analysis_debug_" . session_id() . ".log";
if (!file_exists($debug_log)) {
    touch($debug_log);
    chmod($debug_log, 0664);
}

// Step 2: Handle POST requests in a single switch
if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    switch (true) {
        case isset($_POST['start_session']):
            $name = $_POST['name'] ?? '';
            $surname = $_POST['surname'] ?? '';
            if ($name && $surname) {
                $json_input = json_encode(['action' => 'start_session', 'name' => $name, 'surname' => $surname]);
                $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                    'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
                ]));
                file_put_contents($debug_log, "Login: $json_input\nResponse: " . ($response ?: 'No response') . "\n", FILE_APPEND);
                $data = json_decode($response, true);
                if ($response && isset($data['session_id'])) {
                    $_SESSION['session_id'] = $data['session_id'];
                    $logged_in = true;
                    $results = "Logged in as: " . $_SESSION['session_id'];
                } else {
                    $results = "Login failed: " . ($data['error'] ?? 'Unknown error');
                }
            } else {
                $results = "Name and surname required.";
            }
            break;

        case isset($_POST['logout']):
            $session_id = $_SESSION['session_id'] ?? '';
            if ($session_id) {
                $json_input = json_encode(['action' => 'end_session', 'session_id' => $session_id]);
                @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                    'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
                ]));
                session_unset();
                session_destroy();
                header("Location: analysis_tool.php?logged_out=1");
                exit;
            }
            break;

        case isset($_POST['fetch_proteins']) && $logged_in:
            $protein = $_POST['protein'] ?? '';
            $taxonomy = $_POST['taxonomy'] ?? '';
            $session_id = $_SESSION['session_id'] ?? '';
            if ($protein && $taxonomy && $session_id) {
                $json_input = json_encode(['action' => 'search_protein', 'protein_name' => $protein, 'taxonomy' => $taxonomy, 'session_id' => $session_id]);
                $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                    'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
                ]));
                file_put_contents($debug_log, "Fetch: $json_input\nResponse: " . ($response ?: 'No response') . "\n", FILE_APPEND);
                $data = json_decode($response, true);
                if ($response && !isset($data['error'])) {
                    $results = "<h2>Search Results</h2>";
                    foreach ($data['sequences'] as $seq) {
                        $results .= "<div class='sequence'><strong>{$seq['id']}</strong><br>{$seq['description']}<br><p>Organism: {$seq['organism']}</p><pre>{$seq['sequence']}</pre></div>";
                    }
                    $_SESSION['last_sequence_data'] = json_encode([
                        'protein_name' => $protein,
                        'taxonomy' => $taxonomy,
                        'sequences' => $data['sequences']
                    ]);
                } else {
                    $results = "Error: " . ($data['error'] ?? 'Failed to fetch proteins.');
                }
            } else {
                $results = "Error: Please provide protein and taxonomy.";
            }
            break;

        case isset($_POST['fetch_proteins']) && !$logged_in:
            $results = "Error: Please log in to fetch proteins.";
            break;
    }
}

// Step 3: Fetch user history if logged in
// https://stackoverflow.com/questions/2823913/how-do-i-create-a-user-history 
if ($logged_in) {
    $session_id = $_SESSION['session_id'];
    $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
        'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => json_encode(['action' => 'get_user_history', 'session_id' => $session_id])]
    ]));
    $data_history = json_decode($response, true);
    if ($response && isset($data_history['history'])) {
        $history = $data_history['history'];
    } else {
        $results .= "<p>History fetch failed: " . ($data_history['error'] ?? 'Unknown error') . "</p>";
    }
}

?>

<!DOCTYPE html>
<html>
<head>
    <title>Protein Fetcher - Analysis Tool</title>
    <link rel="stylesheet" href="styles_new.css">
    <style>
        .header-container { display: flex; align-items: center; justify-content: flex-start; gap: 10px; margin-bottom: 20px; }
        .header-container button { padding: 8px 16px; font-size: 14px; }
        .sequence { margin: 10px 0; }
        .history { margin-top: 20px; }
        .history-item { border: 1px solid #ccc; padding: 10px; margin-bottom: 10px; }
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
    <h1>Analysis Tool</h1>

    <?php if (!$logged_in): ?>
        <div id="login">
            <h2>Login</h2>
            <p>Provide your details for an optimized experience.</p>
            <?php if (isset($_GET['logged_out'])): ?>
                <p>Logged out successfully.</p>
            <?php endif; ?>
            <form method="POST">
                <div class="input-group">
                    Name: <input name="name" type="text" placeholder="Enter your name"><br>
                </div>
                <div class="input-group">
                    Surname: <input name="surname" type="text" placeholder="Enter your surname"><br>
                </div>
                <button type="submit" name="start_session">Start</button>
            </form>
            <p><?php echo $results; ?></p>
        </div>
    <?php else: ?>
        <div id="logout">
            <form method="POST">
                <button type="submit" name="logout">Log Out</button>
            </form>
        </div>
        <div id="query">
            <h2>Query Proteins</h2>
            <p>Enter protein and taxonomy to analyze.</p>
            <form method="POST">
                <div class="input-group">
                    Protein: <input name="protein" type="text" placeholder="e.g., insulin"><br>
                </div>
                <div class="input-group">
                    Taxonomy: <input name="taxonomy" type="text" placeholder="e.g., Homo sapiens"><br>
                </div>
                <button type="submit" name="fetch_proteins">Fetch</button>
            </form>
            <form method="GET" action="example_results.php">
                <button type="submit">Try Example (Glucose-6-Phosphatase - Aves)</button>
            </form>
        </div>
        <div id="results">
            <?php echo $results; ?>
            <?php if (isset($data) && !empty($data['sequences'])): ?>
                <form method="POST" action="analyze.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data']); ?>'>
                    <button type="submit">Analyze Conservation</button>
                </form>
                <form method="POST" action="motifscan.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data']); ?>'>
                    <button type="submit">Scan Motifs</button>
                </form>
                <form method="POST" action="advanced.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars($_SESSION['last_sequence_data']); ?>'>
                    <button type="submit">Advanced Analysis</button>
                </form>
            <?php endif; ?>
        </div>
        <?php if (!empty($history)): ?>
            <div class="history">
                <h2>Your Previous Analyses</h2>
                <?php foreach ($history as $item): ?>
                    <div class="history-item">
                        <p><strong>Protein:</strong> <?php echo htmlspecialchars($item['protein_name']); ?></p>
                        <p><strong>Taxonomy:</strong> <?php echo htmlspecialchars($item['taxonomy']); ?></p>
                        <p><strong>Queried:</strong> <?php echo htmlspecialchars($item['query_time']); ?></p>
                        <p><strong>Results:</strong> <?php echo nl2br($item['results']); ?></p>
                    </div>
                <?php endforeach; ?>
            </div>
        <?php endif; ?>
    <?php endif; ?>
</body>
</html>
