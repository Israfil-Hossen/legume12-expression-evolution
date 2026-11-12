# Export the graphical abstract slide with PowerPoint itself (PNG and TIFF at 300 dpi, plus PDF).
# Usage (repository root): powershell -NoProfile -ExecutionPolicy Bypass -File pipeline/11_manuscript/export_graphical_abstract.ps1
$src = (Resolve-Path "manuscript_mpe/figures/GraphicalAbstract.pptx").Path
$dir = Split-Path $src
$pp = New-Object -ComObject PowerPoint.Application
try {
    $deck = $pp.Presentations.Open($src, $true, $false, $false)   # read-only, no window
    $w = [int]($deck.PageSetup.SlideWidth / 72 * 300); $h = [int]($deck.PageSetup.SlideHeight / 72 * 300)
    $deck.Slides(1).Export((Join-Path $dir "GraphicalAbstract.png"), "PNG", $w, $h)
    $deck.Slides(1).Export((Join-Path $dir "GraphicalAbstract.tif"), "TIF", $w, $h)
    $deck.SaveAs((Join-Path $dir "GraphicalAbstract.pdf"), 32)   # ppSaveAsPDF
    $deck.Close()
    "exported ${w} x ${h} px"
} finally { $pp.Quit() }
