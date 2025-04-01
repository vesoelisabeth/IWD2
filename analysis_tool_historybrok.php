<?php
session_start();

// --- Setup and Debugging ---
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);
$debug_log = "/tmp/analysis_debug_" . session_id() . ".log";  // Use /tmp for write access
if (!file_exists($debug_log)) {
    touch($debug_log);
    chmod($debug_log, 0664);
}
file_put_contents($debug_log, "Script started\n", FILE_APPEND);

// --- Session and Variable Initialization ---
$logged_in = isset($_SESSION['session_id']);
$results = '';
$data = null;

// --- Load Sequence Data ---
if (isset($_POST['sequence_data']) && !empty($_POST['sequence_data'])) {
    $data = json_decode($_POST['sequence_data'], true);
    file_put_contents($debug_log, "Loaded sequence data from POST: " . $_POST['sequence_data'] . "\n", FILE_APPEND);
} elseif (isset($_SESSION['last_sequence_data']) && !empty($_SESSION['last_sequence_data'])) {
    $data = json_decode($_SESSION['last_sequence_data'], true);
    file_put_contents($debug_log, "Loaded sequence data from session: " . $_SESSION['last_sequence_data'] . "\n", FILE_APPEND);
}

// --- Handle POST Requests ---
if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    if (isset($_POST['start_session'])) {
        $name = $_POST['name'] ?? '';
        $surname = $_POST['surname'] ?? '';
        if ($name && $surname) {
            $json_input = json_encode(['action' => 'start_session', 'name' => $name, 'surname' => $surname]);
            $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
            ]));
            file_put_contents($debug_log, "Login input: $json_input\nResponse: " . ($response ?: 'No response') . "\n", FILE_APPEND);
            if ($response === false) {
                $results = "Error: Failed to connect to sqlconnect.php.";
            } else {
                $data_response = json_decode($response, true);
                if (isset($data_response['session_id']) && $data_response['session_id']) {
                    $_SESSION['session_id'] = $data_response['session_id'];
                    $logged_in = true;
                    $results = "Session started successfully! Session ID: " . $_SESSION['session_id'];
                } else {
                    $results = "Error: No valid session_id received. Response: " . htmlspecialchars($response);
                }
            }
        } else {
            $results = "Error: Name and surname are required.";
        }
    } elseif (isset($_POST['logout'])) {
        $session_id = $_SESSION['session_id'] ?? '';
        $json_input = json_encode(['action' => 'end_session', 'session_id' => $session_id]);
        $context = stream_context_create([
            'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
        ]);
        @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, $context);
        session_unset();
        session_destroy();
        header("Location: analysis_tool.php?logged_out=1");
        exit;
    } elseif (isset($_POST['fetch_proteins']) && $logged_in) {
        $protein = $_POST['protein'] ?? '';
        $taxonomy = $_POST['taxonomy'] ?? '';
        $session_id = $_SESSION['session_id'] ?? '';
        if ($protein && $taxonomy && $session_id) {
            $json_input = json_encode(['action' => 'search_protein', 'protein_name' => $protein, 'taxonomy' => $taxonomy, 'session_id' => $session_id]);
            $response = @file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
                'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => $json_input]
            ]));
            file_put_contents($debug_log, "Fetch input: $json_input\nResponse: " . ($response ?: 'No response') . "\n", FILE_APPEND);
            if ($response === false) {
                $results = "Error: Failed to connect to sqlconnect.php.";
            } else {
                $data = json_decode($response, true);
                if (!isset($data['error'])) {
                    $results = "<h2>{$data['result']}</h2>";
                    foreach ($data['sequences'] as $seq) {
                        $results .= "<div class='sequence'><strong>{$seq['id']}</strong><br>{$seq['description']}<br><p>Organism: {$seq['organism']}</p><pre>{$seq['sequence']}</pre></div>";
                    }
                    $_SESSION['last_sequence_data'] = json_encode($data);
                } else {
                    $results = "Error: {$data['error']}";
                }
            }
        } else {
            $results = "Error: Please provide protein and taxonomy.";
        }
    } elseif (isset($_POST['fetch_proteins']) && !$logged_in) {
        $results = "Error: Please log in first.";
    }
  }

// --- Fetch history if logged in ---

if ($logged_in) {
    $session_id = $_SESSION['session_id'];
    $response = file_get_contents('https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/sqlconnect.php', false, stream_context_create([
        'http' => ['method' => 'POST', 'header' => 'Content-Type: application/json', 'content' => json_encode(['action' => 'get_user_history', 'session_id' => $session_id])]
    ]));
    $data = json_decode($response, true);
    if (isset($data['history'])) {
        $history = $data['history'];
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
            gap: 10px;
            margin-bottom: 20px;
        }
        .header-container button {
            padding: 8px 16px;
            font-size: 14px;
        }
        .sequence { margin: 10px 0; }
	.history { margin-top: 20px; }
        .history-item { border: 1px solid #ccc; padding: 10px; margin-bottom: 10px; }
    </style>
</head>
<body>
    <!-- Navigation -->
    <nav class="top-nav">
        <a href="indexx.html"><button>Home</button></a>
        <a href="about.php"><button>About</button></a>
        <a href="analysis_tool.php"><button>Analysis Tool</button></a>
        <a href="help.php"><button>Help</button></a>
        <a href="credits.php"><button>Credits</button></a>
    </nav>

    <!-- Page Header -->
    <h1>Analysis Tool</h1>

    <!-- Login Section -->
    <?php if (!$logged_in): ?>
        <div id="login">
            <h2>Login</h2>
            <p>You need to provide your details for optimised user experience.</p>
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
        </div>
    <?php else: ?>
        <!-- Logout Button -->
        <div id="logout">
            <form method="POST">
                <button type="submit" name="logout">Log Out</button>
            </form>
        </div>

        <!-- Protein Query Form -->
        <div id="query">
            <h2>Query Proteins</h2>
            <p>Insert your protein and taxonomy of interest and we will tell you more about it.</p>
            <form method="POST">
                <div class="input-group">
                    Protein: <input name="protein" type="text" placeholder="e.g., insulin"><br>
                </div>
                <div class="input-group">
                    Taxonomy: <input name="taxonomy" type="text" placeholder="e.g., Homo sapiens"><br>
                </div>
                <button type="submit" name="fetch_proteins">Fetch</button>
            </form>
            <!-- Example Button -->
            <form method="GET" action="example_results.php">
                <button type="submit">Try Example (Glucose-6-Phosphatase - Aves)</button>
            </form>
        </div>
        </div>

        <!-- Results and Analysis Options -->
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
                <form method="POST" action="advanced.php">
                    <input type="hidden" name="sequence_data" value='<?php echo htmlspecialchars(json_encode($data)); ?>'>
                    <button type="submit">Advanced Analysis</button>
                </form>
            <?php endif; ?>
        </div>
	<!-- History Section -->
        <?php if (!empty($history)): ?>
            <div class="history">
                <h2>Your Previous Analyses</h2>
                <?php foreach ($history as $item): ?>
                    <div class="history-item">
                        <p><strong>Protein:</strong> <?php echo htmlspecialchars($item['protein_name']); ?></p>
                        <p><strong>Taxonomy:</strong> <?php echo htmlspecialchars($item['taxonomy']); ?></p>
                        <p><strong>Queried:</strong> <?php echo htmlspecialchars($item['query_time']); ?></p>
                        <p><strong>Results:</strong> <?php echo $item['results'] ? htmlspecialchars($item['results']) : 'No results'; ?></p>
                    </div>
	<?php endforeach; ?>
            </div>
        <?php endif; ?>
    <?php endif; ?>
</body>
</html>
