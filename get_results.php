<?php
// Helper file for analysis.html to get the session data 

session_start();
header('Content-Type: application/json');

$alignment = $_SESSION['alignment'] ?? 'No alignment data available';
unset($_SESSION['alignment']); // Clear after use
echo json_encode(['alignment' => $alignment]);
exit;
?>
