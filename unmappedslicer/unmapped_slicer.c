#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <curl/curl.h>

#define BLOCK_SIZE      65536
#define READS_PER_CHUNK 765

curl_off_t content_length = 0;
bool needs_padding = true;

static size_t
write_data (void *ptr, size_t size, size_t nmemb, void *stream)
{
  return fwrite (ptr, size, nmemb, (FILE *)stream);
}

static int
xferinfo (void *p, curl_off_t dltotal, curl_off_t dlnow, curl_off_t ultotal, curl_off_t ulnow)
{
  if (content_length == 0)
    content_length = dltotal;

  if (dlnow > BLOCK_SIZE)
    return 1;

  return 0;
}

static size_t
write_truncated_data (void *ptr, size_t size, size_t nmemb, void *stream)
{
  if (needs_padding == false)
    return fwrite (ptr, size, nmemb, (FILE *)stream);

  uint8_t *data = ptr;
  int32_t byte = 0;

  /*  We look ahead three bytes. */
  for (; byte < (size * nmemb) - 3; byte++)
    {
      /* The BGZF header contains the following byte sequence: 31 139 8 4 0 0 0 0. */
      if (data[byte + 0] == 31 && data[byte + 1] == 139 && data[byte + 2] == 8 && data[byte + 3] == 4)
        {
          printf ("Found start of block at offset %d.\n", byte);
          break;
        }
    }

  int write_bytes = nmemb - byte;
  if (write_bytes > 0)
    {
      fwrite (&(data[byte]), 1, write_bytes, (FILE *)stream);
      needs_padding = false;
    }

  return nmemb;
}


int
main (int argc, char **argv)
{
  CURL *curl_handle;
  CURLcode res;

  FILE *output_file;

  if(argc < 4) {
    printf ("Usage: %s <input-url.bam> <input-url.bam.bai> <out.bam>\n", argv[0]);
    return 1;
  }

  curl_global_init (CURL_GLOBAL_ALL);

  char *input_bam  = argv[1];
  char *input_bai  = argv[2];
  char *output_bam = argv[3];

  size_t output_bai_len = strlen (argv[3]) + 5;
  char output_bai[output_bai_len];
  if (snprintf (output_bai, output_bai_len, "%s.bai", output_bam) < 0)
    exit (1);

  /* 
   * Download the header.
   * ------------------------------------------------------------------------ */

  curl_handle = curl_easy_init ();

  output_file = fopen (output_bam, "wb");
  if(! output_file)
    {
      fprintf (stderr, "Cannot open '%s'\n", output_bam);
      exit (1);
    }

  /* FIXME: Remove these. */
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0L);
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0L);

  curl_easy_setopt (curl_handle, CURLOPT_URL, input_bam);
  curl_easy_setopt (curl_handle, CURLOPT_XFERINFOFUNCTION, xferinfo);
  curl_easy_setopt (curl_handle, CURLOPT_NOPROGRESS, 0L);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEFUNCTION, write_data);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEDATA, output_file);

  /* The following assumes that reading one block will fetch the
   * entire header.  This is a 64kb download. */
  curl_easy_perform (curl_handle);

  long end_position = ftell (output_file);
  fseek (output_file, 2L, SEEK_SET);
  uint8_t *buffer = calloc (1, end_position - 2);
  long byte = 0;
  fread (buffer, 1, end_position - 2, output_file);
  for (; byte < end_position - 5; byte++)
    {
      /* The BGZF header contains the following byte sequence: 31 139 8 4 0 0 0 0. */
      if (buffer[byte + 0] == 31 && buffer[byte + 1] == 139 && buffer[byte + 2] == 8 && buffer[byte + 3] == 4)
        {
          fprintf (stderr, "Found first data block at offset %d.\n", byte);
          break;
        }
    }

  /* Truncate the file at a new BGZF block. */
  if (byte < end_position - 5)
    ftruncate (fileno (output_file), end_position - byte);

  fclose (output_file);
  curl_easy_cleanup (curl_handle);

  fprintf (stderr, "Downloaded BAM header.\n");

  /* 
   * Download the BAI.
   * ------------------------------------------------------------------------ */

  curl_handle = curl_easy_init ();

  output_file = fopen (output_bai, "wb");
  if(! output_file)
    {
      fprintf (stderr, "Cannot open '%s'\n", output_bai);
      exit (1);
    }

  /* FIXME: Remove these. */
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0L);
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0L);

  curl_easy_setopt (curl_handle, CURLOPT_XFERINFOFUNCTION, NULL);
  curl_easy_setopt (curl_handle, CURLOPT_NOPROGRESS, 1L);

  curl_easy_setopt (curl_handle, CURLOPT_URL, input_bai);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEFUNCTION, write_data);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEDATA, output_file);
  curl_easy_perform (curl_handle);
  res = curl_easy_perform (curl_handle);
  if (res != CURLE_OK)
    fprintf (stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror (res));

  fclose (output_file);
  curl_easy_cleanup (curl_handle);
  fprintf (stderr, "Downloaded BAI file.\n");

  /*
   * Read the BAI + header.
   * ------------------------------------------------------------------------ */

  BGZF* bam_stream = NULL;
  hts_idx_t *bam_idx = NULL;
  bam_hdr_t *bam_hdr = NULL;

  bam_hdr = bam_hdr_init ();
  bam_stream = bgzf_open (output_bam, "r");
  bam_hdr = bam_hdr_read (bam_stream);
  bam_idx = bam_index_load (output_bam);

  if (bam_idx == NULL)
    {
      fprintf (stderr, "Cannot read the BAM index.\n");
      exit (1);
    }

  uint64_t unmapped_reads = 0;
  uint64_t mapped_reads = 0;
  uint64_t total_reads = 0;
  int32_t chromosome = 0;
  for (; chromosome < bam_hdr->n_targets; chromosome++)
    {
      hts_idx_get_stat (bam_idx, chromosome, &mapped_reads, &unmapped_reads);
      printf ("Chromosome %s: %lu mapped, %lu unmapped.\n",
              bam_hdr->target_name[chromosome],
              mapped_reads,
              unmapped_reads);

      total_reads += mapped_reads;
      total_reads += unmapped_reads;

      if (unmapped_reads > 0)
        printf ("Chromosome %d has %lu unmapped reads\n", chromosome, unmapped_reads);
    }

  bgzf_close (bam_stream);

  /* One chunk has a maximum size of 2^16 (that's 64 * 1024). */
  unmapped_reads = hts_idx_get_n_no_coor (bam_idx);
  uint64_t read_after = total_reads / READS_PER_CHUNK * 64 * 1024;

  uint64_t offset_position = content_length / total_reads - BLOCK_SIZE * 4;
  double percentage = (float)unmapped_reads / (float)total_reads;
  printf ("Unmapped reads: %lu (%f%%)\n", unmapped_reads, percentage * 100);
  printf ("It's almost safe to read after %lu bytes.\n", read_after);

  /* 
   * Download the unmapped reads.
   * ------------------------------------------------------------------------ */

  curl_handle = curl_easy_init ();

  output_file = fopen (output_bam, "ab");
  if(! output_file)
    {
      fprintf (stderr, "Cannot open '%s'\n", output_bam);
      exit (1);
    }

  /* The last 28 bytes are the BGZF null block.  Position the file pointer such that
   * we overwrite this block. */
  fseek (output_file, -28, SEEK_END);
  long size = ftell (output_file);
  ftruncate (fileno (output_file), size);

  char *range = calloc (1, 255);
  snprintf (range, 255, "%lu-%" CURL_FORMAT_CURL_OFF_T, read_after, read_after + content_length);
  fprintf (stderr, "Going to download %" CURL_FORMAT_CURL_OFF_T " megabytes\n",
           (content_length - read_after) / 1024 / 1024);

  /* FIXME: Remove these. */
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0L);
  curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0L);
  //curl_easy_setopt (curl_handle, CURLOPT_VERBOSE, 1L);

  curl_easy_setopt (curl_handle, CURLOPT_URL, input_bam);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEFUNCTION, write_truncated_data);
  curl_easy_setopt (curl_handle, CURLOPT_WRITEDATA, output_file);
  curl_easy_setopt (curl_handle, CURLOPT_RANGE, range);

  res = curl_easy_perform (curl_handle);
  if (res != CURLE_OK)
    fprintf (stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror (res));

  /* Always close the file properly. */
  char end[] = { 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
                 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
                 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  fwrite (end, 1, 28, output_file);
  
  free (range);
  fclose (output_file);
  curl_easy_cleanup (curl_handle);
  fprintf (stderr, "Downloaded unmapped reads.\n");

  curl_global_cleanup ();
  return 0;
}
