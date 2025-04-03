<?php
session_start();

$file = $_GET['file'] ?? '';

if (empty($file)) {
    http_response_code(400);
    die("Invalid request: No file specified.");
}

$path = urldecode($file); // Decode URL-encoded path

if (file_exists($path) && (strpos($path, '/home/s2015320/public_html/project2/tmp/motifscan_') === 0 || 
                           strpos($path, '/home/s2015320/public_html/project2/tmp/advanced_') === 0 || 
                           strpos($path, '/home/s2015320/public_html/project2/tmp/conservation_') === 0 ||
			   strpos($path, '/home/s2015320/public_html/project2/tmp/alignment_') === 0)) {
    header('Content-Type: text/plain');
    header('Content-Disposition: attachment; filename="' . basename($path) . '"');
    readfile($path);

    // Cleanup

    if (is_file($path)) {
        unlink($path); // Remove only the alignment file, no chmod or dir cleanup
    }
} else {
    http_response_code(404);
    echo "File not found: " . htmlspecialchars($path);
}
exit;
?>
