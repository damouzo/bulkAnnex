# palettes.R — shared colour helpers for bulkAnnex pipeline scripts.
#
# ba_distinct_pal() returns a *deterministic* set of well-separated colours for a
# set of categorical levels, using only base R (grDevices) — no new packages.
# Levels are sorted internally so that the same set of levels always maps to the
# same colours, independent of call order. This keeps the pipeline PNGs and the
# dashboard in sync (a given sample/condition always gets the same colour).
#
# Strategy: small sets use curated qualitative palettes; larger sets fall back to
# a golden-angle hue wheel that maximises distance between consecutive colours.
# Note: no palette can make more than ~10-12 categories truly distinguishable.
ba_distinct_pal <- function(levels) {
    levels <- as.character(levels)
    levels <- sort(unique(levels[!is.na(levels) & nzchar(levels)]))
    n <- length(levels)
    cols <- if (n == 0) {
        character(0)
    } else if (n <= 8) {
        grDevices::hcl.colors(n, palette = "Dark 3")
    } else if (n <= 12) {
        grDevices::hcl.colors(n, palette = "Set 3")
    } else {
        grDevices::hcl(h = (seq_len(n) * 137.508) %% 360, c = 55, l = 65)
    }
    stats::setNames(cols, levels)
}
