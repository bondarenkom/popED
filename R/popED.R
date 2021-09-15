#' @title Population Distributed Equally
#' 
#' @author Maksym Bondarenko <mb4@soton.ac.uk>,
#'         David Kerr <dk2n16@soton.ac.uk>,  
#'         Alessandro Sorichetta <as1v13@soton.ac.uk> and
#'         Tom McKeen <t.l.mckeen@soton.ac.uk>   
#'         
#' @details This function produces gridded population equally distributed
#' 
#'  
#' @usage
#' popED(mastergrid, px_area, pop, cores = 4, output_dir=tempdir(), 
#' prefix="A1", nt=2, check_result=TRUE, verbose = TRUE)
#' 
#' @param mastergrid Raster tat contains the unique area IDs as their value. 
#' @param px_area Raster containing the pixel area. 
#' @param pop Character vector containing the name of the file from which the 
#'        unique area ID and corresponding population values are to be read 
#'        from. The file should contain two columns comma-separated with the 
#'        value of administrative ID and population without columns names. 
#'        If it does not contain an absolute path, the file name is relative to
#'        the current working directory.
#' @param cores Integer vector containing an integer. Indicates the number of 
#'        cores to use in parallel when executing the function. If set to 0 
#'        \code{(max_number_of_cores - 1)}  will be used based on as many
#'        processors as the hardware and RAM allow. Default is \code{cores} = 0.
#' @param output_dir Character vector containing the path to the directory for 
#'        writing output files. Default is the temp directory.
#' @param prefix Character to be added to the output results. Default is "A1"
#' @param nt Numeric. Increase number of blocks sugesting for 
#'        processing raster file. Default is \code{nt} = 2
#' @param check_result Logical vector indicating whether the results will be 
#'        compared with input data. Default is \code{check_result} = TRUE.
#' @param verbose Logical vector indicating whether to print 
#'        intermediate output from the function to the console, which might be 
#'        helpful for model debugging. Default is \code{verbose} = TRUE.
#' @importFrom raster nlayers raster
#' @importFrom utils write.csv
#' @rdname popED
#' @return Raster* object of gridded population equally distributed.
#' @export
#' @examples
#' \dontrun{
#' 
#' library("popED")
#' 
#' res <- popED(mastergrid="input_mastergrid.tif", 
#'              px_area="input_px_area.tif", 
#'              pop="pop_table.csv",
#'              cores=4,
#'              output_dir="/user/output" 
#'              ) 
#'  
#' # Plot populataion raster 
#' plot(res) 
#' 
#' }
popED <- function(mastergrid,
                  px_area,
                  pop,
                  cores = 4,
                  output_dir=tempdir(),
                  prefix="A1",
                  nt=2,
                  check_result=TRUE,
                  verbose = TRUE){
  
  if(!file.exists( mastergrid )) {
    
    msg <- paste0("Error ::  'mastergrid' file does not exist")
    stop(msg)
    
  } 
  
  if(!file.exists( px_area )) {
    
    msg <- paste0("Error :: 'px_area' file does not exist")
    stop(msg)
    
  } 
  
  if(!file.exists( pop )) {
    
    msg <- paste0("Error ::  'pop' file does not exist")
    stop(msg)
    
  } 
  
  if (is.null(output_dir) | !dir.exists(output_dir)) {
    
    msg <- paste0("Error :: Output directory ",output_dir," does not exsit. Please choose a different directory.")
    stop(msg)
    
  }  
  
  
  #f_b <- basename(mastergrid)
  mastergrid_base_name = substr(basename(mastergrid), 1, nchar(basename(mastergrid)) - 4)
  px_area_base_name = substr(basename(px_area), 1, nchar(basename(px_area)) - 4)
  
  mastergrid.rst <- raster(mastergrid)
  px_area.rst <- raster(px_area)
  
  df <- utils::read.csv(pop, stringsAsFactors=FALSE, header = FALSE)
  colnames(df) <-  c("ADMINID", "ADMINPOP") 
  
  blocks_size <- get_blocks_size(mastergrid.rst, cores, n=nt)      
  
  npoc_blocks <- ifelse(blocks_size$n < cores, blocks_size$n, cores)
  
  #
  log_info("MSG", paste0("Rasterizing population table using ", mastergrid_base_name ), 
           verbose=verbose)
  
  rst.pop.census <- rasterize_parallel(mastergrid.rst, 
                                       df, 
                                       cores=npoc_blocks, 
                                       blocks=blocks_size,
                                       filename=file.path(output_dir, paste0("rasterized_cencus_",prefix,".tif")), 
                                       overwrite=TRUE, 
                                       silent=F)
  
  
  log_info("MSG", paste0("Calculating zonal stats  ", px_area_base_name ), 
           verbose=verbose)
  
  PixelArea.zonal.stats.sum <- calculate_zs_parallel(px_area.rst,
                                                     mastergrid.rst,
                                                     fun='sum', 
                                                     cores=npoc_blocks, 
                                                     blocks=blocks_size, 
                                                     na.rm=TRUE, 
                                                     silent=TRUE)  
  
  
  if (verbose){
    cat("\n")  
  }
  
  log_info("MSG", paste0("Rasterizing pixel area  "), 
           verbose=verbose)
  
  rst.PixelArea.zonal.stats <- rasterize_parallel(mastergrid.rst, 
                                                  PixelArea.zonal.stats.sum, 
                                                  cores=npoc_blocks, 
                                                  blocks=blocks_size,
                                                  filename=file.path(output_dir, paste0("rasterized_px_area_",prefix,".tif")), 
                                                  overwrite=TRUE, 
                                                  silent=F)
  
  
  log_info("MSG", paste0("Calculating pop distributed equally  "), 
           verbose=verbose)
  
  pop_pixel <- (px_area.rst * rst.pop.census)/rst.PixelArea.zonal.stats
  
  
  log_info("MSG", paste0("Writing results to file  ", file.path(output_dir, paste0("pop_dist_equally_",prefix,".tif"))), 
           verbose=verbose)
  
  writeRaster(pop_pixel, 
              filename=file.path(output_dir, paste0("pop_dist_equally_",prefix,".tif")),
              format="GTiff", 
              datatype='FLT4S', 
              overwrite=TRUE, 
              options=c("COMPRESS=LZW"),
              NAflag=-99999)
  
  
  if (check_result){
    
    zonal.stats.sum <- calculate_zs_parallel(pop_pixel,
                                             mastergrid.rst,
                                             fun='sum', 
                                             cores=npoc_blocks, 
                                             blocks=blocks_size, 
                                             na.rm=TRUE, 
                                             silent=TRUE) 
    
    
    colnames(zonal.stats.sum) <- c("ADMINID", "PPP_FINAL_RES")
    output_stats.sorted <-  zonal.stats.sum[sort.list(zonal.stats.sum[,1]), ] 
    output_stats.sorted <-  output_stats.sorted[output_stats.sorted[,1] != 0, ] 
    
    df_f <- merge( as.data.frame(df), 
                   as.data.frame(output_stats.sorted), 
                   by="ADMINID", 
                   sort=FALSE) 
    
    
    df_f <- cbind( df_f, abs(df_f$ADMINPOP - df_f$PPP_FINAL_RES)) 
    
    colnames(df_f) <- c("ADMINID", "ADMINPOP", "PPP", "DIFFERENCE")
    
    file.path.csv <- file.path(output_dir, 
                               paste0("check_result_",prefix,".csv"))
    
    write.csv( as.data.frame(df_f), file = file.path.csv, row.names=FALSE )
    
    
    log_info("MSG", paste0("Checking  "), 
             verbose=verbose)
    cat("------------------------------------------------\n")
    
    log_info("MSG", paste0("Accuracy/error is  ", mean(df_f[,4])), 
             verbose=verbose)
    cat("------------------------------------------------\n")
    
    
  }
  
  
}
