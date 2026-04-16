process SHINY_DASHBOARD {
    tag "dashboard"
    label 'process_single'

    container params.container

    publishDir "${params.outdir}", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    val sync  // sync signal: number of upstream versions.yml collected (no file staging)

    output:
    path 'dashboard/',   emit: dashboard_dir
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p dashboard/modules dashboard/www

    cp    ${projectDir}/dashboard/app.R               dashboard/
    cp    ${projectDir}/dashboard/global.R            dashboard/
    cp    ${projectDir}/dashboard/ui.R                dashboard/
    cp    ${projectDir}/dashboard/server.R            dashboard/
    cp -r ${projectDir}/dashboard/modules/.           dashboard/modules/
    cp -r ${projectDir}/dashboard/www/.               dashboard/www/
    cp    ${projectDir}/dashboard/launch_dashboard.sh     dashboard/ 2>/dev/null || true
    cp    ${projectDir}/dashboard/launch_dashboard_hpc.sh dashboard/ 2>/dev/null || true

    Rscript -e "writeLines(c('SHINY_DASHBOARD:', paste0('    R: ', R.version[['version.string']]), paste0('    shiny: ', as.character(packageVersion('shiny')))), 'versions.yml')"
    """
}
