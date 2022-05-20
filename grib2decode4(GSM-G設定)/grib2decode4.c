/*----------------------------------------------------------------------
  2007.08.22		 sample_decoder.c

  GRIB2 ファイル サンプルデコード処理プログラム

	このプログラムの全部又は一部を利用してもかまいませんが、利用したこ
  とによって利用者が被った直接的又は間接的ないかなる損害についても、
  気象庁は一切責任を負いません。
	また、プログラムに関する個別の対応は行いかねますので、ご容赦願います。

  利用方法

  ANSI 準拠の C コンパイラでコンパイルして下さい。
  GRIB2ファイルのファイル名を引数に与えて実行すると、
  ファイルの内容が表示されます。

  更新履歴
  2018.02.02 LFMフォーマットに対応
  2018.09.18 潮位フォーマットに対応、切り出し機能追加
  2019.02.15 抽出対象のパラメータカテゴリー、パラメータ番号指定機能追加
  2019.03.01 切り出しを緯度経度指定ではなくXY番号指定に変更
  2020.04.23 GSMガイダンス/MSMガイダンス/GSM_GPV/週間アンサンブル フォーマットに対応
----------------------------------------------------------------------*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <string.h>
# include <time.h>

/* 1,2,4,8 : signed integer   (x byte)
		 u : unsigned integer (1 byte)
		 S : unsigned integer (2 byte)
		 C : character		(4 byte)
		 R : float			(4 byte) */

char *sectionFormat[9] = {
  "C2uu8",			 /* Section 0 */
  "4u22uuu2uuuuuuu",   /* Section 1 */
  "",				  /* Section 2 */
  "4uu4uuS",		 /* + Section 3 */
  "4u2S",			/* + Section 4 */
  "4u4S",			/* + Section 5 */
  "4uu",			 /* + Section 6 */
  "4u",			  /* + Section 7 */
  "C"				  /* Section 8 */
};

typedef struct { int secno, templat; char *format; } TemplateFormat;

TemplateFormat  templateFormat[] = {
  { 3, 0, "uu4u4u4444444u4444u" },
  { 3,20, "uu4u4u44444u4444uu" },
  { 4, 0, "uuuuu2uu4u14u14" },
  { 4, 1, "uuuuu2uu4u14u14uuu" },
  { 4, 8, "uuuuu2uu4u14u142uuuuuu4uuu4u4" },
  { 4, 9, "uuuuu2uu4u14u14uuu14142uuuuuu4uuu4u4" },
  { 4,11, "uuuuu2uu4u14u14uuu2uuuuuu4uuu4u4" },
  { 5, 0, "R22uu" },
  { 0, 0, "" }
};

static int isLittleEndian;

# define FIT_BYTE_ORDER(pointer,size) if(isLittleEndian)swabN(pointer,size,1)

/* reverse byte order */
static void
swabN( void *buf, int size, int nn )
{
	char *ba, *bb, *buf2 = buf;
	while( nn-- ) {
		bb = ( ba = buf2 ) + size -1;
		do {
			char a;
			a   = *ba;
			*ba = *bb;
			*bb =  a;
		} while( ++ba < --bb );
		buf2 += size;
	}
}

int read_section_0( FILE *fp, void **sec_buffer )
{
	char *bufr;
	*sec_buffer = realloc( *sec_buffer, 16 );
	bufr = *sec_buffer;
	if( fread( bufr, 1, 16, fp ) != 16 ||
		strncmp( bufr, "GRIB", 4 ) != 0 ) {
		fprintf( stderr, "Really GRIB file ?\n" );
		exit(1);
	}
	return 0;
}

int read_section_X( FILE *fp, void **sec_buffer )
{
	int  length;
	unsigned char *bufp;

	fread( &length, 4, 1, fp );
	if( strncmp( (char *)(&length), "7777", 4 ) == 0 ) {
		*sec_buffer = realloc( *sec_buffer, 4 );
		strncpy( *sec_buffer, "7777", 4 );
		return 8;
	}
	FIT_BYTE_ORDER( &length, 4 );

	*sec_buffer = realloc( *sec_buffer, length );
	bufp = *sec_buffer;
	if( fread( bufp + 4, 1, length - 4, fp ) != length - 4 ) {
		fprintf( stderr, "Unexpected EOF\n" );
		exit(1);
	}
	FIT_BYTE_ORDER( &length, 4 );
	memcpy( bufp, &length, 4 );

	return (int)bufp[4];   /* section No. */
}

void decode_buf( char *sec_buffer, int *index, char *format, double **dvpp )
{
	unsigned char buffer[8];
			  int  size, ii, jj, missing;
		   double  *dvp, dd;

	dvp = *dvpp;

	for( ii = *index ; *format ; format++ ) {
		switch( *format ) {
			case '1' :
			case '2' :
			case '4' :
			case '8' : size = *format - '0'; break;
			case 'u' : size = 1; break;
			case 'S' : size = 2; break;
			case 'R' :
			case 'C' : size = 4; break;
			default :
				fprintf( stderr, "Internal Error !\n" );
				exit(9);
		}
		memcpy( buffer, sec_buffer+ii, size );

		missing = 1;
		for( jj = 0 ; jj < size ; jj++ )
			missing &= ( buffer[jj] == 0xFF );

		switch( *format ) {
			case 'u' : dd = *buffer;					  break;
			case 'S' : dd = (buffer[0]<<8) + buffer[1];   break;

			case '1' :
			case '2' :
			case '4' :
			case '8' : dd = buffer[0] & 0x7F;
					   for( jj = 1 ; jj < size ; jj++ )
						   dd = dd * 256.0 + buffer[jj];
					   if( *buffer & 0x80 ) dd = -dd;	 break;

			case 'R' : FIT_BYTE_ORDER( buffer, size );
					   dd = *(float *)buffer;			 break;

			case 'C' : memcpy( &dd, buffer, size );	   break;
		}
		//if( size == 1 ) printf( "   %4d	 : ", ii+1 );
		//else			printf( "%4d .. %-4d: ", ii+1, ii+size );

		//	 if( missing )		printf( "missing \n" );
		//else if( *format == 'R' ) printf( "%f\n", dd );
		//else if( *format == 'C' ) printf( "'%.4s'\n", &dd );
		//else					  printf( "%.0f\n", dd );

		*dvp++ = dd;
		ii += size;
	}
	*index = ii;
	*dvpp = dvp;
}

void decode_section( int secno, char *sec_buffer, double **double_values )
{
	int index, length;
	double *dvp;
	char *format;
	TemplateFormat *tfp;

	switch( secno ) {
	case 0 : length =  5; break;
	case 1 : length = 15; break;
	case 7 : length =  2; break;
	case 8 : length =  1; break;
	default :
		memcpy( &length, sec_buffer, 4 );
		FIT_BYTE_ORDER( &length, 4 );
	}
	dvp = realloc( *double_values, sizeof(double) * length );
	*double_values  = dvp;

	//printf( "== SECTION %d == \n", secno );

	format = sectionFormat[secno];
	index = 0;

	decode_buf( sec_buffer, &index, format, &dvp );

	if( 3 <= secno && secno <= 5 ) {
		int templat = (int)dvp[-1];
		for( tfp = templateFormat ; tfp->secno ; tfp++ ) {
			if( tfp->secno == secno &&
				tfp->templat == templat ) break;
		}
		if( ! tfp->secno ) {
			fprintf( stderr, "No Information about Template %d.%d\n",
					secno, templat );
			return;
		}
		decode_buf( sec_buffer, &index, tfp->format, &dvp );
	}
}

void unpack_data( const char *sec_7, const double *values_5, float **out )
{
	double  rr, ee, dd, pow_2_e, pow_10_d;
	float  *op;
	int	 num, nbit, ii;
	unsigned int uw, mask;
	const unsigned char *uc;

	uc = (const unsigned char *)sec_7 + 5;

	num   = values_5[2];
	rr	= values_5[4];
	ee	= values_5[5];
	dd	= values_5[6];
	nbit  = values_5[7];
	pow_2_e  = pow( 2.0, ee);
	pow_10_d = pow(10.0, dd);

	*out = realloc( *out, sizeof(**out) * num );
	op = *out;

	mask = 1;
	for( ii = 0 ; ii < nbit - 1 ; ii++ )
		mask = mask << 1 | 1;

	for( ii = 0 ; ii < num ; ii++ ) {
		memcpy( &uw, &uc[(nbit*ii)/8], 4 );
		FIT_BYTE_ORDER( &uw, 4 );
		uw = ( uw>>(32-nbit-(nbit*ii)%8)) & mask;
		*op++ = ( rr + pow_2_e * uw ) / pow_10_d;
	}
}


void unpack_bitmap( const char *sec_6, const double *values_3,  unsigned char **out )
{
	unsigned char  *op;
	int	 num, nbit, ii;
	unsigned int uw, mask;
	const unsigned char *uc;

	uc = (const unsigned char *)sec_6 + 6;

	num   = values_3[3];
	nbit = 1;

	*out = realloc( *out, sizeof(**out) * num );
	op = *out;

	mask = 1;
	for( ii = 0 ; ii < nbit - 1 ; ii++ )
		mask = mask << 1 | 1;

	for( ii = 0 ; ii < num ; ii++ ) {
		memcpy( &uw, &uc[(nbit*ii)/8], 1 );

		uw = ( uw>>(8-nbit-(nbit*ii)%8)) & mask;
		*op++ = uw;
	}
}



void save_file( char *outPath,char *outName, const double *values_3, const double *values_4, const double *values_5,
		const unsigned char *bitmap_data,
		const float *unpacked_data,
		int cut_mode, int *cut_lat, int *cut_lon,
		int yy, int mm, int dd, int hh, int mi, int ss,
		int prodfld, int datakind)
{
	FILE *fpout;
	char   fileName[256];
	char *elem;
	int ii;
	int xx;
	int ix, pos;
	int num;
	int paramc, paramn, fixedk, fornum, termnum, fortm, enskind, ensno, tempno;
	double oa;
	char *paramc_nm, *paramn_nm, *fixedk_nm;
	int  sta_y, end_y, sta_x, end_x;
	int  xsize,ysize, s_lat, s_lon, e_lat, e_lon, dx, dy;
	int  ic, cut_x[11], cut_y[11];
	double  cut_ff[11];

	tempno = (int)values_4[3];   //プロダクト定義テンプレート番号
	fornum = (int)values_4[12];  //予報時間
	enskind = -1;
	if(tempno == 9){
		termnum = (int)values_4[37]; //統計処理した期間
	}else if(tempno == 11){
		enskind = (int)values_4[19]; //19アンサンブル予報の種類
		ensno = (int)values_4[20];   //20摂動番号
		termnum = (int)values_4[33]; //統計処理した期間
	}else{
		termnum = (int)values_4[30]; //統計処理した期間
	}
	fortm = termnum + fornum;        //_予測時間（予報時間＋統計処理した期間）

	fprintf( stderr, "template no 4.%d ,fortime=%2d, ruika=%2d  \n", tempno, fornum, termnum );

	if(enskind >= 0)
		sprintf( fileName, "%s%s_M%d%02d_%03d.csv", outPath, outName, enskind, ensno, fortm ); //_Mメンバー_予測時間
	else
		sprintf( fileName, "%s%s_LEN%02d_%02d.csv", outPath, outName, termnum, fortm ); //_LEN統計処理した期間_予測時間
	//else if(datakind > 0 && datakind < 100)
	//	sprintf( fileName, "%s%s_LEN%02d_%02d.csv", outPath, outName, termnum, fornum ); //_LEN統計処理した期間_予報時間
	//else
	//	sprintf( fileName, "%s%s_%02d.csv", outPath, outName, fornum ); //_予報時間

	fprintf( stdout, "output '%s'\n", fileName );

	if( ( fpout = fopen( fileName, "wb" ) ) == NULL ) {
		fprintf( stderr, "file '%s' open error! \n", fileName );
		exit(1);
	}

	num = (int)values_5[2];		/* 全資料点の数 */
	paramc = (int)values_4[4];	/* パラメータカテゴリー */
	paramn = (int)values_4[5];	/* パラメータ番号 */
	fixedk = (int)values_4[13];	/* 第一固定面の種類 */
	oa = -2; 		/* 範囲外の場合の値 */

	if(prodfld == 0){	//0:気象プロダクト
		//upd_s 20190215
		//elem = "TotalPrecipitation";
		if(paramc == 1 && paramn == 8)
			elem = "TotalPrecipitation"; //総降水量
		else
			elem = "";
		//upd_e 20190215

		//パラメータカテゴリー
		if(paramc == 0)
			paramc_nm = "温度";
		else if(paramc == 1)
			paramc_nm = "湿度";
		else if(paramc == 2)
			paramc_nm = "運動量";
		else if(paramc == 3)
			paramc_nm = "質量";
		else if(paramc == 4)
			paramc_nm = "短波放射";
		else
			paramc_nm = "";
		//パラメータ番号
		if(paramc == 0 && paramn == 0)
			paramn_nm = "温度(K)";
		else if(paramc == 1 && paramn == 1)
			paramn_nm = "相対湿度(％)";
		else if(paramc == 1 && paramn == 8)  //add 20200423
			paramn_nm = "総降水量";
		else if(paramc == 1 && paramn == 52) //add 20200423
			paramn_nm = "降水強度の合計";
		else if(paramc == 2 && paramn == 2)
			paramn_nm = "風のｕ成分(m/s)";
		else if(paramc == 2 && paramn == 3)
			paramn_nm = "風のｖ成分(m/s)";
		else if(paramc == 2 && paramn == 8)
			paramn_nm = "鉛直速度（気圧）(Pa/s)";
		else if(paramc == 3 && paramn == 0)
			paramn_nm = "気圧(Pa/s)";
		else if(paramc == 3 && paramn == 1)
			paramn_nm = "海面更正気圧(Pa)";
		else if(paramc == 3 && paramn == 5)
			paramn_nm = "ジオポテンシャル高度(gpm)";
		else if(paramc == 4 && paramn == 7)
			paramn_nm = "下向き短波放射フラックス(Wm)";
		else
			paramn_nm = "";
	}
	else if(prodfld == 10){	//10:海洋プロダクト
		oa = -999;
		elem = "";
		//パラメータカテゴリー
		if(paramc == 3)
			paramc_nm = "海表面の特性";
		else
			paramc_nm = "";
		//パラメータ番号
		if(paramn == 1)
			paramn_nm = "平均海面からの偏差";
		else if(paramn == 200)
			paramn_nm = "天文潮位";
		else
			paramn_nm = "";
	}else{
		elem = "";
		paramc_nm = "";
		paramn_nm = "";
	}
	//第一固定面の種類
	if(fixedk == 1)
		fixedk_nm = "地面";
	else if(fixedk == 101)
		fixedk_nm = "平均海面";
	else if(fixedk == 103)
		fixedk_nm = "地上からの特定高度面";
	else if(fixedk == 100)
		fixedk_nm = "等圧面(Pa)";
	else
		fixedk_nm = "";

	xsize = (int)values_3[14];	//X格子点数
	ysize = (int)values_3[15];	//Y格子点数
	s_lat = (int)values_3[18];	//最初の格子点の緯度
	s_lon = (int)values_3[19];	//最初の格子点の経度
	e_lat =  (int)values_3[21];	//最後の格子点の緯度
	e_lon =  (int)values_3[22];	//最後の格子点の経度
	dx = (int)values_3[23];		//i方向の増分
	dy = (int)values_3[24];		//j方向の増分

	sta_y = 0;
	sta_x = 0;
	end_y = ysize - 1;
	end_x = xsize - 1;

	for( ic = 0 ; ic < sizeof(cut_y); ic++ ) {
		cut_ff[ic] = oa;
		cut_y[ic] = -1;
		cut_x[ic] = -1;
	}

	if(cut_mode == 0){
		//del_s 20190301
		//if(cut_lat[0] <= s_lat && cut_lon[0] >= s_lon && cut_lat[1] >= e_lat && cut_lon[1] <= e_lon) {
		//	//切り取り範囲指定
		//	sta_y = (ysize - 1) - (int)((cut_lat[0] - e_lat) / (double)dy);
		//	sta_x = (int)((cut_lon[0] - s_lon) / (double)dx);

		//	end_y = (ysize - 1) - (int)((cut_lat[1] - e_lat) / (double)dy);
		//	end_x = (int)((cut_lon[1] - s_lon) / (double)dx);

		//	e_lat = e_lat + (dy * ((ysize - 1) - end_y));
		//	e_lon = s_lon + (dx * end_x);
		//	s_lat = s_lat - (dy * sta_y);
		//	s_lon = s_lon + (dx * sta_x);
		//}
		//del_e 20190301
		//add_s 20190301
		if(cut_lat[0] > 0)
			sta_y = cut_lat[0] - 1;
		if(cut_lon[0] > 0)
			sta_x = cut_lon[0] - 1;
		if(cut_lat[1] <=  ysize)
			end_y = cut_lat[1] - 1;
		if(cut_lon[1] <=  xsize)
			end_x = cut_lon[1] - 1;
		//add_e 20190301

		//add_s 2020.04.23
		fprintf( stdout, "area %d %d %d %d\n", s_lat, s_lon,  e_lat, e_lon );
		fprintf( stdout, "cut area %d %d %d %d\n", sta_y, sta_x,  end_y, end_x );

		e_lat = e_lat + (dy * ((ysize - 1) - end_y));
		e_lon = s_lon + (dx * end_x);
		s_lat = s_lat - (dy * sta_y);
		s_lon = s_lon + (dx * sta_x);

		fprintf( stdout, "cut area %d %d %d %d\n", s_lat, s_lon,  e_lat, e_lon );
		//add_e 2020.04.23
	}else if(cut_mode > 0){
		for( ic = 0 ; ic < cut_mode; ic++ ) {
			cut_y[ic] = cut_lat[ic];
			cut_x[ic] = cut_lon[ic];
		}
	}

	/* ヘッダ出力 */
	fprintf( fpout, "資料の参照時刻          	 , %d/%02d/%02d %02d:%02d:%02d, %s\n", yy, mm, dd, hh, mi, ss, elem);
	fprintf( fpout, "最初の格子点の緯度      	 , %d\n", s_lat );
	fprintf( fpout, "最初の格子点の経度      	 , %d\n", s_lon );
	fprintf( fpout, "最後の格子点の緯度      	 , %d\n", e_lat );
	fprintf( fpout, "最後の格子点の経度      	 , %d\n", e_lon );
	fprintf( fpout, "i方向の増分             	 , %d\n", dx );
	fprintf( fpout, "j方向の増分             	 , %d\n", dy );
	fprintf( fpout, "パラメータカテゴリー    	 , %d, %s\n", paramc, paramc_nm );
	fprintf( fpout, "パラメータ番号          	 , %d, %s\n", paramn, paramn_nm );
	fprintf( fpout, "第一固定面の種類        	 , %d, %s\n", fixedk, fixedk_nm );
	fprintf( fpout, "第一固定面の尺度因子    	 , %d\n", (int)values_4[14] );
	fprintf( fpout, "第一固定面の尺度付きの値	 , %d\n", (int)values_4[15] );
	fprintf( fpout, "全資料点の数            	 , %d\n", num );
	fprintf( fpout, "参照値(R)               	 , %f\n", (double)values_5[4] );
	fprintf( fpout, "二進尺度因子(E)         	 , %f\n", (double)values_5[5] );
	fprintf( fpout, "十進尺度因子(D)         	 , %f\n", (double)values_5[6] );
	fprintf( fpout, "範囲外とする値          	 , %.2f,", (double)oa);	//add 20180202
	if(cut_mode == 0)
		fprintf( fpout, " %d,%d,%d,%d, //開始x,y,終了x,y(1,1を基点)", sta_x + 1, sta_y + 1, end_x + 1, end_y + 1);	//add 20180202 upd 20190301
	if(cut_mode > 0){
		for( ic = 0 ; ic < cut_mode; ic++ ) {
			fprintf( fpout, " %d:%d,", cut_lat[ic], cut_lon[ic]);
		}
	}
	fprintf( fpout, "\n");

	ix = 0;
	xx = 0;
	double ff=0;
	for( ii = 0 ; ii < ysize; ii++ ) {
		for (xx=0; xx<xsize; xx++) {

			ff = oa;
			pos=(ii*xsize+xx);
			if( bitmap_data == NULL || bitmap_data[pos] == 1 ){
				ff = unpacked_data[ix++];
			}

			if(cut_mode > 0){
				for( ic = 0; ic < cut_mode; ic++ ) {
					if(cut_y[ic] == ii && cut_x[ic] == xx){
						cut_ff[ic] = ff;
					}
				}
			}else{

				if(ii >= sta_y && ii <= end_y && xx >= sta_x && xx <= end_x){
					if(xx < end_x){
						fprintf( fpout, "%.2f,", ff );
					}else{
						// x格子点数まできたら改行
						fprintf( fpout, "%.2f\n", ff );
					}
				}
			}

		}
	}
	if(cut_mode > 0){
		for( ic = 0; ic < cut_mode; ic++ ) {
			if(ic < (cut_mode - 1)){
				fprintf( fpout, "%.2f,", cut_ff[ic]  );
			}else{
				fprintf( fpout, "%.2f\n", cut_ff[ic]  );
			}
		}
	}


	fclose( fpout );
}


int main( int argc, char *argv[] )
{
	FILE *fpin;
	char *fileName;
	int	 secno, ii;
	int	 cnt;
	void   *sec_buffer = NULL;
	float  *unpacked_data = NULL;
	unsigned char *bitmap_data = NULL;
	double *sec_double_value[9];

	int  ic, cut_mode;
	int  cut_lat[11], cut_lon[11];
	int  datakind,target_flg,prodfld,paramc, paramn, fixedk, cremdl;
	int tgt_kind, tgt_kind1, tgt_kind2, tempno4;	//add 20190215


	// check system byte order
	isLittleEndian = 1;
	isLittleEndian = *(char *)(&isLittleEndian);

	for( ii = 0 ; ii < sizeof(sec_double_value)/sizeof(*sec_double_value) ;  )
		sec_double_value[ii++] = NULL;
	if( argc == 1 ) {
		 fprintf( stderr, "\n\n usage: %s 'grib2 file name'\n\n", argv[0] );
		 exit(1);
	}

	char *outPath = "";
	char *outName = "";
	if (argc>2){
		outPath=argv[2];
		outName=argv[3];
	}

	tgt_kind = 0;	//add 20190215
	cut_mode = -1;	//0:XY範囲切り出し, 1～:抽出XYの数
	for(ic=0;ic<sizeof(cut_lat);ic++){
		cut_lat[ic] = 0;
		cut_lon[ic] = 0;
	}
	if (argc>4){
		cut_mode = atoi(argv[4]);
	}
	if (cut_mode == 0){
		//切り出し開始X、Y、終了X、Y
		if(argc>8){
			//del_s 20190301
			//cut_lat[0] = atof(argv[5]) * pow(10.0, 6);
			//cut_lon[0] = atof(argv[6]) * pow(10.0, 6);
			//cut_lat[1] = atof(argv[7]) * pow(10.0, 6);
			//cut_lon[1] = atof(argv[8]) * pow(10.0, 6);
			//del_e 20190301
			//add_s 20190301
			cut_lon[0] = atoi(argv[5]) ; //x
			cut_lat[0] = atoi(argv[6]);  //y
			cut_lon[1] = atoi(argv[7]) ; //x
			cut_lat[1] = atoi(argv[8]);  //y
			//add_e 20190301
		}else{
			fprintf( stderr, "\n\n usage: %d 'cut out x y range'\n\n", cut_mode ); //upd 20190301
			exit(1);
		}
		if(argc>9){	//add 20190215
			tgt_kind=atoi(argv[9]);
		}
	}else if(cut_mode > 0){
		//指定したxy地点を抽出
		if(cut_mode > 10)
			cut_mode = 10;
		if(argc>(5+cut_mode)){
			ii=5;
			for(ic=0;ic<cut_mode;ic++){
				cut_lon[ic] = atoi(argv[ii++]) ; //x
				cut_lat[ic] = atoi(argv[ii++]);  //y
			}
		}else{
			fprintf( stderr, "\n\n usage: %d 'pic up x y'\n\n", cut_mode );
			exit(1);
		}
	}
	else{
		if(argc>5){	//add 20190215
			tgt_kind=atoi(argv[5]);
		}
	}


	fileName = argv[1];
	if( ( fpin = fopen( fileName, "rb" ) ) == NULL ) {
		fprintf( stderr, "grib2 file '%s' open error! \n", fileName );
		exit(1);
	}

	datakind = 0;	//0:LFM/MSM, 100:tide
	if(strstr(fileName, "_TID_")!=NULL ||  strstr(fileName, "_SGM_")!=NULL)
		datakind = 100;
	else if(strstr(fileName, "_MSM_GUID_")!=NULL)
		datakind = 2;
	else if(strstr(fileName, "_GSM_GUID_")!=NULL)
		datakind = 3;
	else if(strstr(fileName, "_GSM_GPV_")!=NULL)
		datakind = 4;
	else if(strstr(fileName, "_EPSG_GPV_")!=NULL)
		datakind = 21;    //週間アンサンブル

	if(tgt_kind==0 && (datakind == 2 || datakind == 3)){
		tgt_kind=152;
	}else if(tgt_kind==0){
		tgt_kind=18;
	}

	//add_s 20190215
	if(datakind < 100){
		if(tgt_kind < 100){
			tgt_kind1=(int)(tgt_kind / 10);
			tgt_kind2=(tgt_kind % 10);
		}else{
			tgt_kind1=(int)(tgt_kind / 100);
			tgt_kind2=(tgt_kind % 100);
		}
	}
	//add_e 20190215

	secno = read_section_0( fpin, &sec_buffer );
	decode_section( secno, sec_buffer, &sec_double_value[secno] );

	cnt = 0;
	while( secno = read_section_X( fpin, &sec_buffer ) ) {
		decode_section( secno, sec_buffer, &sec_double_value[secno] );
		if( secno == 8 ) break;
		if( secno == 6  && (int)sec_double_value[6][2] == 0) {
			unpack_bitmap( sec_buffer, sec_double_value[3], &bitmap_data );
		}
		if( secno == 7 ) {
			target_flg = 0;
			prodfld = 0;
			tempno4 = (int)sec_double_value[4][3];		//product template no
			paramc = (int)sec_double_value[4][4];		//parameter category
			paramn = (int)sec_double_value[4][5];		//parameter number
			fixedk = (int)sec_double_value[4][13];		//type of fixed surface
			cremdl = (int)sec_double_value[4][7];		//model
			if(datakind < 100){      //LFM/MSM
				prodfld = 0;	//0:気象プロダクト
				//1:湿度, 8:総降水量
				//if(paramc == 1 && paramn == 8)	//del 20190215
				if(paramc == tgt_kind1 && paramn == tgt_kind2)	//upd 20190215
					target_flg = 1;
				if(tempno4 == 9)
					target_flg = 0;
			}
			else if(datakind == 100) 	{    //tide
				prodfld = 10;	//10:海洋プロダクト
				//1:地面, 3:海表面の特性,  1:平均海面からの偏差/200:天文潮位,   225:高潮モデル/226:天文潮位モデル
				if(fixedk == 1 && paramc == 3 && (paramn == 1 || paramn == 200) && (cremdl == 225 || cremdl == 226))
					target_flg = 1;
			}

			if (target_flg==1){
				unpack_data( sec_buffer, sec_double_value[5], &unpacked_data );

				save_file ( outPath,outName, sec_double_value[3], sec_double_value[4], sec_double_value[5],
								bitmap_data, unpacked_data,
								cut_mode, cut_lat, cut_lon,
								 (int)sec_double_value[1][7],	/* 年 */
								 (int)sec_double_value[1][8],	/* 月 */
								 (int)sec_double_value[1][9],	/* 日 */
								 (int)sec_double_value[1][10],	/* 時 */
								 (int)sec_double_value[1][11],	/* 分 */
								 (int)sec_double_value[1][12],	/* 秒 */
								 prodfld,datakind
				);
			}
		}
	}
	fclose(fpin);
	return 0;
}

