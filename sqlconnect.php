<?php
ob_start();
ini_set('display_errors', 1);
error_reporting(E_ALL);

$localhost = "127.0.0.1";
$username = "s2015320";
$password = "Qk3UxizAsaJO8ld_w6xrb1xdG2HRXo";
$dbname = "s2015320_website";

header('Content-Type: application/json');

try {
    $conn = new PDO("mysql:host=$localhost;dbname=$dbname", $username, $password);
    $conn->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);

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
    }

    $input = json_decode(file_get_contents('php://input'), true);
    $action = $input['action'] ?? '';

    if ($action === 'init_db') {
        init_db($conn);
        ob_end_clean();
        echo json_encode(["result" => "Database initialised"]);
        exit;
    }

    if ($action === 'start_session') {
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
    }

    // Terminate the session option
    if ($action === 'end_session') {
        $session_id = $input['session_id'] ?? '';
        if ($session_id) {
            $stmt = $conn->prepare("DELETE FROM users WHERE session_id = ?");
            $stmt->execute([$session_id]);
            // CASCADE will remove related protein_queries and protein_results
            ob_end_clean();
            echo json_encode(["result" => "Session ended"]);
        } else {
            ob_end_clean();
            echo json_encode(["error" => "No session_id provided"]);
        }
        exit;
    }

    // Searching for the protein 
    if ($action === 'search_protein') {
        $protein_name = $input['protein_name'] ?? '';
        $taxonomy = $input['taxonomy'] ?? '';
        $session_id = $input['session_id'] ?? ''; 
        
        if (!$protein_name || !$taxonomy || !$session_id) {
            ob_end_clean();
            echo json_encode(["error" => "Missing required fields"]);
            exit;
        }

        // Store the query. Step 1.4
        try {
            $stmt = $conn->prepare("INSERT INTO protein_queries (session_id, protein_name, taxonomy, query_time) VALUES (?, ?, ?, NOW())");
            $stmt->execute([$session_id, $protein_name, $taxonomy]);
        } catch (PDOException $e) {
            ob_end_clean();
            die(json_encode(["error" => "Database error: " . $e->getMessage()]));
        }

        // Connecting to the BioPython script
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
            die(json_encode([
                "error" => "Python script failed",
                "command" => $command,
                "output" => $result,
                "return_code" => $return_var
            ]));
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

        $stmt = $conn->prepare("
            INSERT INTO protein_results (session_id, protein_name, taxonomy, protein_id, description, sequence, organism)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ");
        foreach ($data['sequences'] as $seq) {
            $stmt->execute([$session_id, $protein_name, $taxonomy, $seq['id'], $seq['description'], $seq['sequence'], $seq['organism']]);
        }
        
        ob_end_clean();
        echo json_encode($data);
        exit;
    }

    ob_end_clean();
    echo json_encode(["error" => "Invalid action"]);
    $conn = null;

    // w3school https://www.w3schools.com/php/php_mysql_connect.asp
    // BPSM Biopython lecture
    // SQL lectures from BPSM and IWD2
    // https://stackoverflow.com/questions/52535136/python-return-json-to-php-and-use-json-decode
} catch(PDOException $e) {
    ob_end_clean();
    die(json_encode(["error" => "Connection failed: " . $e->getMessage()]));
}
?>
