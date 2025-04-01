<?php
ob_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

// Database connection settings
$localhost = "127.0.0.1";
$username = "s2015320";
$password = "Qk3UxizAsaJO8ld_w6xrb1xdG2HRXo";
$dbname = "s2015320_website";

header('Content-Type: application/json');

try {
    $conn = new PDO("mysql:host=$localhost;dbname=$dbname", $username, $password);
    $conn->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);

    // Table creation function
    function init_db($conn) {
        $conn->exec("CREATE TABLE IF NOT EXISTS users (
            id INT AUTO_INCREMENT PRIMARY KEY,
            name VARCHAR(50),
            surname VARCHAR(50),
            session_id VARCHAR(50) UNIQUE,
            login_time DATETIME
        )");
        $conn->exec("CREATE TABLE IF NOT EXISTS protein_queries (
            id INT AUTO_INCREMENT PRIMARY KEY,
            session_id VARCHAR(50),
            protein_name VARCHAR(100),
            taxonomy VARCHAR(100),
            query_time DATETIME,
            FOREIGN KEY (session_id) REFERENCES users(session_id) ON DELETE CASCADE
        )");
        $conn->exec("CREATE TABLE IF NOT EXISTS protein_results (
            id INT AUTO_INCREMENT PRIMARY KEY,
            session_id VARCHAR(50),
            protein_name VARCHAR(100),
            taxonomy VARCHAR(100),
            protein_id VARCHAR(50),
            description TEXT,
            sequence TEXT,
            organism VARCHAR(100),
            FOREIGN KEY (session_id) REFERENCES users(session_id) ON DELETE CASCADE
        )");
        $conn->exec("CREATE TABLE IF NOT EXISTS jobs (
            job_id VARCHAR(50) PRIMARY KEY,
            session_id VARCHAR(50),
            files TEXT,
            timestamp DATETIME
        )");
    }

    // Parse incoming JSON request
    $input = json_decode(file_get_contents('php://input'), true);
    $action = $input['action'] ?? '';

    // Initialize database tables
    if ($action === 'init_db') {
        init_db($conn);
        ob_end_clean();
        echo json_encode(["result" => "Database initialised"]);
        exit;
    // Handle user login
    } elseif ($action === 'start_session') {
        $name = $input['name'] ?? '';
        $surname = $input['surname'] ?? '';
        if (!$name || !$surname) {
            ob_end_clean();
            echo json_encode(["error" => "Name and surname required"]);
            exit;
        }
        $session_id = $name . "_" . $surname . "_" . time();
        $stmt = $conn->prepare("INSERT INTO users (name, surname, session_id, login_time) VALUES (?, ?, ?, NOW())");
        $stmt->execute([$name, $surname, $session_id]);
        ob_end_clean();
        echo json_encode(["session_id" => $session_id]);
        exit;
    // Handle user logout
    } elseif ($action === 'end_session') {
        $session_id = $input['session_id'] ?? '';
        if ($session_id) {
            $stmt = $conn->prepare("DELETE FROM users WHERE session_id = ?");
            $stmt->execute([$session_id]);
            ob_end_clean();
            echo json_encode(["result" => "Session ended"]);
        } else {
            ob_end_clean();
            echo json_encode(["error" => "No session_id provided"]);
        }
        exit;
    // Search proteins using Biopython connection
    } elseif ($action === 'search_protein') {
        $protein_name = $input['protein_name'] ?? '';
        $taxonomy = $input['taxonomy'] ?? '';
        $session_id = $input['session_id'] ?? '';
        if (!$protein_name || !$taxonomy || !$session_id) {
            ob_end_clean();
            echo json_encode(["error" => "Missing required fields"]);
            exit;
        }
        $stmt = $conn->prepare("INSERT INTO protein_queries (session_id, protein_name, taxonomy, query_time) VALUES (?, ?, ?, NOW())");
        $stmt->execute([$session_id, $protein_name, $taxonomy]);
        $json_input = json_encode(["protein_name" => $protein_name, "taxonomy" => $taxonomy]);
        $script_path = "/home/s2015320/public_html/project2/biopython_connect.py";
        $python_path = "/usr/bin/python3";
        $temp_file = tempnam(sys_get_temp_dir(), 'bio_');
        file_put_contents($temp_file, $json_input);
        $command = "$python_path $script_path < $temp_file 2>&1";
        exec($command, $output, $return_var);
        $result = implode("\n", $output);
        unlink($temp_file);
        if ($return_var !== 0 || empty($result)) {
            ob_end_clean();
            die(json_encode(["error" => "Python script failed", "command" => $command, "output" => $result, "return_code" => $return_var]));
        }
        $result = preg_replace('/^Content-type: application\/json\s*/i', '', $result);
        $data = json_decode($result, true);
        if ($data === null) {
            ob_end_clean();
            die(json_encode(["error" => "Invalid JSON from Python script: " . $result]));
        }
        if (isset($data['error'])) {
            ob_end_clean();
            echo json_encode($data);
            exit;
        }
        $stmt = $conn->prepare("INSERT INTO protein_results (session_id, protein_name, taxonomy, protein_id, description, sequence, organism) VALUES (?, ?, ?, ?, ?, ?, ?)");
        foreach ($data['sequences'] as $seq) {
            $stmt->execute([$session_id, $protein_name, $taxonomy, $seq['id'], $seq['description'], $seq['sequence'], $seq['organism']]);
        }
        ob_end_clean();
        echo json_encode($data);
        exit;
    // Store job results
    } elseif ($action === 'store_job') {
        $job_id = $input['job_id'] ?? '';
        $session_id = $input['session_id'] ?? '';
        $files = $input['files'] ?? [];
        if (!$job_id || !$session_id || empty($files)) {
            ob_end_clean();
            echo json_encode(["error" => "Missing required fields"]);
            exit;
        }
        $stmt = $conn->prepare("INSERT INTO jobs (job_id, session_id, files, timestamp) VALUES (?, ?, ?, NOW()) ON DUPLICATE KEY UPDATE files = ?, timestamp = NOW())");
        $stmt->execute([$job_id, $session_id, json_encode($files), json_encode($files)]);
        ob_end_clean();
        echo json_encode(["success" => true, "job_id" => $job_id]);
        exit;
    // Retrieve job details
    } elseif ($action === 'get_job') {
        $job_id = $input['job_id'] ?? '';
        if (!$job_id) {
            ob_end_clean();
            echo json_encode(["error" => "No job_id provided"]);
            exit;
        }
        $stmt = $conn->prepare("SELECT * FROM jobs WHERE job_id = ?");
        $stmt->execute([$job_id]);
        $result = $stmt->fetch(PDO::FETCH_ASSOC);
        if ($result) {
            $result['files'] = json_decode($result['files'], true);
            ob_end_clean();
            echo json_encode($result);
        } else {
            ob_end_clean();
            echo json_encode(["error" => "Job not found"]);
        }
        exit;
    // Fetch user history
    } elseif ($action === 'get_user_history') {  // New action
        $session_id = $input['session_id'] ?? '';
        if (!$session_id) {
            ob_end_clean();
            echo json_encode(["error" => "No session_id provided"]);
            exit;
        }
        $stmt = $conn->prepare("
            SELECT pq.protein_name, pq.taxonomy, pq.query_time, 
                   GROUP_CONCAT(CONCAT(pr.protein_id, ': ', pr.description) SEPARATOR '<br>') as results
            FROM protein_queries pq
            LEFT JOIN protein_results pr ON pq.session_id = pr.session_id 
                AND pq.protein_name = pr.protein_name 
                AND pq.taxonomy = pr.taxonomy
            WHERE pq.session_id = ?
            GROUP BY pq.protein_name, pq.taxonomy, pq.query_time
            ORDER BY pq.query_time DESC
        ");
        $stmt->execute([$session_id]);
        $history = $stmt->fetchAll(PDO::FETCH_ASSOC);
        ob_end_clean();
        echo json_encode(["history" => $history]);
        exit;
    // Handle invalid actions
    } else {
        ob_end_clean();
        echo json_encode(["error" => "Invalid action"]);
        exit;
    }

    // References:
    // w3school https://www.w3schools.com/php/php_mysql_connect.asp
    // BPSM Biopython lecture
    // SQL lectures from BPSM and IWD2
    // https://stackoverflow.com/questions/52535136/python-return-json-to-php-and-use-json-decode
} catch(PDOException $e) {
    ob_end_clean();
    die(json_encode(["error" => "Connection failed: " . $e->getMessage()]));
}
?>
