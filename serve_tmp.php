<?php
$file = $_GET['file'] ?? '';
$path = "/home/s2015320/public_html/project2/tmp/" . basename($file);

if (file_exists($path) && strpos($path, '/home/s2015320/public_html/project2/tmp/') === 0) {
    if (pathinfo($path, PATHINFO_EXTENSION) === 'png') {
        header('Content-Type: image/png');
    } else {
        header('Content-Type: text/plain');
        header('Content-Disposition: attachment; filename="' . basename($path) . '"');
    }
    readfile($path);
//    unlink($path);
} else {
    http_response_code(404);
    echo "File not found: " . htmlspecialchars($path);
}
?>
