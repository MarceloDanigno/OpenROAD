/**************************************************************************
***
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2007 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
***  Original Affiliation:   UCLA, Computer Science Department,
***                          Los Angeles, CA 90095-1596 USA
***
***  Permission is hereby granted, free of charge, to any person obtaining
***  a copy of this software and associated documentation files (the
***  "Software"), to deal in the Software without restriction, including
***  without limitation
***  the rights to use, copy, modify, merge, publish, distribute, sublicense,
***  and/or sell copies of the Software, and to permit persons to whom the
***  Software is furnished to do so, subject to the following conditions:
***
***  The above copyright notice and this permission notice shall be included
***  in all copies or substantial portions of the Software.
***
*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***
***************************************************************************/

//  Created : 10/09/97, Mike Oliver, VLSI CAD ABKGROUP UCLA

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <string.h>
#include "abkrand.h"

static unsigned seedbuf1279[] = {0xd19b0078, 0x80c91c2a, 0xad234a7d, 0xfe4d556f, 0x9a506e9a, 0x5c306c36, 0xf9b0a447, 0xd0a341df, 0x58bc02a7, 0xc2e76ca2, 0x7047bd4a, 0xf4777fe3, 0x270fe60a, 0xdd62ab57, 0xca904318, 0x978971ab, 0x8a690fef, 0x19908c3b, 0x0fc40762, 0xe125ac90, 0xea2cca93, 0xce2250ac, 0x50ff718d, 0xdb5fae9a, 0x3779fb76, 0xd12307cd, 0xea0b6b2f, 0x53d1bdef, 0xa1eae734, 0x0ef12415, 0x8b3fc239, 0xa078bc47, 0x31dc5d99, 0x529e50b0, 0x977572d0, 0x4b0d8f4d, 0x06d0c208, 0x77939570, 0xd4ef7558, 0x4d4063a9, 0xeafb52b8, 0xc9b2f444, 0xf64bfada, 0xf887fcf5, 0xdb6213a4, 0xd395f95a, 0xfceb1c04, 0xb4455a17, 0x03eecf75, 0x3bdfb037, 0xbe49be14, 0x7543d5b2, 0x64785ea2, 0x02aeae03, 0xf3f2e755, 0x76a4e39e, 0xd8ee97f5, 0xbc5c27e8, 0x591b5722, 0x088582d8, 0x89402829, 0xfd3760e0, 0x29142e20, 0xea491041, 0xbca70723, 0x24096026, 0x1379e71b, 0x78609fb5, 0xb0fe3782, 0xd565c9a1, 0x6ab1bd37, 0x19faaaaa, 0xb75857a4, 0x094fe514, 0xba3ca65e, 0x7501a4a8, 0xaa2a067d, 0x960970e9, 0x0a5b9728, 0xceb9a751,
                                 0x793f7159, 0x2dc55afd, 0x58c02d24, 0xe10f856f, 0x37276c87, 0xfe1453d1, 0x60c8f6c8, 0x5ed58663, 0x998d47c7, 0xcada3295, 0x72a1d309, 0xb33e939e, 0x1ceff6ef, 0x30e46553, 0xcafafadd, 0x73e01abc, 0x44d168a1, 0xc0b1c7c5, 0xe407a525, 0xbfaf899a, 0xc64ccb3e, 0x82a526c6, 0xfcd18eaf, 0x111182c5, 0x3ab267b8, 0x44221a12, 0x073e848f, 0x744a4e97, 0x0267ae89, 0x946c0e6c, 0x87ef2d0b, 0x0f222a21, 0xf21028fa, 0x38f4f6c8, 0xf7d9059e, 0x13f3014a, 0xad80ff94, 0xbdcfb9b6, 0x29ab29c2, 0xbf28379f, 0x3a440d0a, 0xca4fa56b, 0xfdd70b33, 0x06a29f82, 0x289220d5, 0x25ea2f50, 0x839bcbf3, 0xbc64f12b, 0x137ff636, 0x20148af9, 0xd6ed13c0, 0x466786ac, 0xcd4cda11, 0xf76eaa10, 0x181a2f6c, 0xa7fcd59a, 0xc3e7cdfb, 0x4b58ac87, 0x8bec3841, 0xe7193a1e, 0x328b3da4, 0x20794a4b, 0x4f26d4ae, 0x0ae38e17, 0x53c84584, 0x4fc6cf52, 0x11dbf4cc, 0x466c7dfe, 0xba0f5ce8, 0xbaf0c09b, 0x34f90a20, 0x619413d9, 0xa6766f45, 0x80236b8a, 0x06fcc5c3, 0x911ff7b2, 0xc9099243, 0xdfbfa8fe, 0xc39c24fc, 0xf2fc5706,
                                 0x0335716c, 0x6fe8eb42, 0x3a3a7d43, 0x7493489d, 0xaaf978ec, 0x563a0331, 0x4b73a5c6, 0xe3df4243, 0x7d1199a2, 0xdd7a3830, 0x5afd20b8, 0x8251c773, 0xb5b1aef1, 0xab1fc224, 0x0e13c97f, 0xb1cf834b, 0x2d048cba, 0x0ce3ec77, 0x4aed344e, 0x5fd165c9, 0xbfeb21b2, 0x85116a4f, 0xc17f7de6, 0x75472d95, 0x63a59b3c, 0x19af6664, 0x6dad2465, 0xbd6cd4f7, 0x333ed161, 0x837c7d38, 0x3ac0325e, 0x89e526ae, 0xdab2ca76, 0x7886fcf5, 0x78b82511, 0x508277d3, 0xb7e640e9, 0x68f5e1a3, 0xcb742fe7, 0x03b1b08d, 0xdbf0aa2c, 0x9734defe, 0x1ba95efd, 0xaa7db64a, 0x369e4170, 0x4ac80ace, 0x68b02c1e, 0x7acaa622, 0xeacb14d9, 0x75555098, 0x1ce352f4, 0xcd6ac2f0, 0x78516b1b, 0xe7d9c1da, 0x7d33b38d, 0x1ca133e3, 0xa3660f56, 0xeaa5f64c, 0xdcdcd4d4, 0xa9c366ec, 0x1e0c1db5, 0xf03b72c2, 0xe7955b99, 0xa1b7ad7a, 0xffcdba94, 0x923b91c7, 0x8e5cf1d0, 0xc427e76c, 0x750cb046, 0xf7a74acf, 0xdaf22cac, 0x1f49d7b5, 0xb1d8f741, 0x71fa097d, 0xc1175fcc, 0x5a094180, 0xb8929037, 0x27a48554, 0x26532b14, 0x61adcbb4,
                                 0x97df02fe, 0x412e5468, 0xdb39fecf, 0xdf76c58c, 0xf5789a23, 0x1a840077, 0xd2492878, 0x2b9ab93b, 0x70bef221, 0x31690dfa, 0x67c0a668, 0x17856031, 0x3e6c2ce6, 0xdc525f79, 0xbc9fdbea, 0x44b728a9, 0xe8475b9d, 0xed47d5bd, 0x6bcb5cd3, 0x99ffa691, 0xdb402936, 0xa8007fbf, 0x1bf4b177, 0x400c2620, 0x91b99bb8, 0x9eb94b61, 0x25539563, 0x46e5b7d5, 0xe2e493b4, 0xb18a6c43, 0x4fd446ba, 0x6d321cb9, 0x9d0b70dc, 0x8be102b8, 0x598277f9, 0x7d86e355, 0x9da6c914, 0x59b3b690, 0xb482925d, 0x315b2c37, 0x1e5ea94f, 0xb6c4f865, 0x24d34702, 0x5fb448b2, 0xeeaffecd, 0x23940746, 0x0abe27c9, 0x9acf213d, 0xeb6a62ca, 0xc17c350f, 0xf8561b33, 0x4510edd8, 0x644ec4fe, 0xaa990298, 0xb19f0e65, 0x946f7e9a, 0xb78360f8, 0x13ded0ca, 0xcaa097c2, 0xad00618f, 0xe53baf96, 0xc19f5a11, 0x094e0d5b, 0xd0a12612, 0x30e9b058, 0x40acf15b, 0x0f2ee8b7, 0x90d927dc, 0xc467fe9f, 0x12c671b1, 0x5358fea6, 0x7854fdc3, 0x91b5480e, 0x3529ee6e, 0xaf6afbdb, 0x05fe50c1, 0x6ddfebe0, 0x9b968470, 0xcec56aeb, 0x18eaae66,
                                 0xb06d5a03, 0xf2444ed9, 0x67cd0320, 0xda06e02a, 0x51b23053, 0xaaa064bf, 0x478505e2, 0x7c6a1410, 0x418841a7, 0x064cee34, 0xa5539ceb, 0x3330fd69, 0x242be16d, 0x089a8e61, 0x807016ca, 0x885e5e57, 0x7956e5b1, 0x2d96be71, 0x9eecc954, 0xc1a986b4, 0x416febf5, 0x2090464b, 0x6a4d83f4, 0x0fed69a6, 0x15a9310a, 0xfa2ff5a0, 0xdc0b33da, 0x2d8bc612, 0x07b51f14, 0x3973ae03, 0xaed73f14, 0x94f18d70, 0xc7c415a5, 0xde612d87, 0x6b7a1659, 0xe014eea0, 0xbe9e24a6, 0x43f10682, 0xd29de2ab, 0x2800ce77, 0x08f84e17, 0x9ce179f8, 0x8271e78e, 0x6ccbd3ef, 0x7c2e2a67, 0xb63e99cb, 0x1bb8d845, 0x6946f7cc, 0xaa3ab13c, 0x09193a17, 0x86963b20, 0x0613d2e9, 0x369d30cd, 0x4487d809, 0x5ee436be, 0x6e711f8a, 0xecc12315, 0x86579d6e, 0x060587e6, 0x7f82dee6, 0x37a0809e, 0x5ddfbc5a, 0x971c4156, 0x9a27e9c3, 0x789a4021, 0x58502032, 0xddc7dbfa, 0xef4b5d1d, 0xc563342a, 0x8c086982, 0x3c8612ec, 0xde20aeb8, 0x297508e9, 0x55297e30, 0xcea68b6a, 0xc6678b7a, 0x7144b8c3, 0xd7f21e51, 0x1d7c1c4e, 0x2da7eb90,
                                 0xaec4db87, 0x86cad972, 0x920252fe, 0x214af5c6, 0xbb2e3752, 0xb6128086, 0x7b3cf055, 0xacde2d1b, 0x81906482, 0x6101874f, 0x26d47d14, 0x074f2cc2, 0x115ce833, 0x9e363ea1, 0x6f9f1e8a, 0x2f423d34, 0x6e94fdc9, 0xe520aaec, 0x5dc4e65c, 0xcc74fd11, 0x25932bcf, 0xe4cef30b, 0x52ddc765, 0x955c8146, 0x5c44c003, 0x68c76424, 0xcd17403b, 0xee700230, 0x03358c69, 0xf29d8447, 0xf4e7dfc7, 0xd27e9f8d, 0x12002e68, 0xbb480753, 0x1461520c, 0xbcea3734, 0x46179ce0, 0x9b68b5c0, 0xb33f60e8, 0x1f5a9c43, 0x7e8d9379, 0xac649ebf, 0x62ab13dd, 0xf040aa24, 0xf3ce3ca9, 0xc2d7e5ff, 0xcf415d8a, 0xb1ec0c38, 0x8e015622, 0xbd6875e4, 0xc4bcd131, 0xef6209bf, 0x655d8748, 0x8f425701, 0xb2cc46dc, 0x595825a9, 0x728e0bc2, 0xda9c5477, 0xf8d9297d, 0xc0792c44, 0x911b96f1, 0x71a487ad, 0xa2fc2fc7, 0xc0bd0a4a, 0xe74cb2eb, 0x8b037e28, 0x2fc58c4d, 0x42a82286, 0x5523a961, 0xce113a4b, 0xdf5f9c9a, 0x02b84fdc, 0x2ba055aa, 0xa93bba14, 0x05efd677, 0x22b264a3, 0x5184ee79, 0x2ca1b100, 0xdfb8e04a, 0x95c5a5e7,
                                 0xee32eaf1, 0xc4dc0e56, 0xaad5bf71, 0x9305b4ca, 0x8a48535c, 0x6a256c06, 0x819b9612, 0xc00d3750, 0x7db35bdc, 0x21fb91c4, 0xd4991ad5, 0xcca64adf, 0x57dde2da, 0xf3340e25, 0xd95c4467, 0x98b8a596, 0x94f7d18d, 0x65ae0367, 0xf09d339d, 0xdaa50d28, 0x2ca4a0c0, 0x35eec8a1, 0x252c7f8a, 0x7c5637ae, 0x22932314, 0x467aa38b, 0x414bcd2f, 0xfe49d15b, 0x60652d76, 0x002f7252, 0x8e887372, 0xa683e22b, 0x20815d9d, 0xb567a893, 0x560f228a, 0xe7283894, 0x81347dca, 0xdf6db091, 0x6d657a71, 0x38932bc3, 0x3753d9ce, 0x0d29808b, 0xaff38318, 0xce740449, 0x65817ff1, 0x5aaf50a6, 0x1adc9d01, 0x99bbd3e3, 0x5f8f4159, 0xbe389a2b, 0x4731d93f, 0xcb87e052, 0x2ccd7e31, 0x4a42c859, 0x8e3c6d76, 0xf836840e, 0x0fde27ab, 0x7b1b9ec2, 0xe9d80b40, 0xb0cc56e8, 0x3052494e, 0xc9c23c6c, 0x34be30ca, 0xc3f4b55c, 0xf1ef75ab, 0x9392fc00, 0x2e0cab55, 0xac70cc58, 0xfc934cb4, 0x064246a5, 0x59501884, 0xd6ca4619, 0xe63b50ce, 0x6c7e657f, 0x5dae50d9, 0x09c947f6, 0x6d679cd2, 0xac3a1a70, 0x7ca0948e, 0xbe1e9ca7,
                                 0x9a81b2b5, 0xdfa9387e, 0x09fff47f, 0x0aafe5da, 0xd70818c1, 0x498254c3, 0x6a7cf352, 0xec423967, 0x5b4fe9d4, 0x172efbe6, 0xcad7c50c, 0xee6e8056, 0x368d7e1a, 0x4d7c6474, 0x7b1f1231, 0x2c1a1cd6, 0xd01bba69, 0x6b2b7a45, 0x13fdf925, 0xc2d0940d, 0xd143ed7b, 0xd2ff2176, 0xee6d828a, 0x385a4e99, 0x9fbdb9ce, 0xa59421e2, 0xd213588b, 0xeb9f93a3, 0xdb82afaa, 0xc8b3e41b, 0xc8f82855, 0x7af3cdc4, 0xeba14683, 0x2e3614f6, 0xf2233312, 0xf02517f8, 0xb1473141, 0xa7cfd85a, 0x02ba2d5d, 0x9e6c31f1, 0xcf967fb1, 0xe4506974, 0xf62fffd2, 0xc9fdcae5, 0x0ada4f9d, 0x4ac9dc33, 0x8dffe025, 0x9c3079e2, 0x1014eb72, 0x65df01d2, 0x31d9c5ab, 0xff907a82, 0x06e6b21a, 0xc49285ee, 0xc38b796e, 0xf7996ad4, 0x94ef8278, 0xbeade129, 0x80a67e66, 0x20f4bd22, 0x778b734a, 0xe9abe8dc, 0x7abee6bd, 0xa95d5183, 0xfa373573, 0x5db5db3c, 0x1b1375a5, 0xf240485a, 0x1378a011, 0x13d5d253, 0x1052367b, 0x9b098363, 0xc0681bb1, 0x608ab174, 0x7c1a20d6, 0xac9bf2cb, 0x2696e2e8, 0x8689a8f2, 0xdb21df3b, 0x1c52df13,
                                 0x7814e96a, 0x424d977b, 0x1c18a6cc, 0x019f8ee6, 0x079e18d8, 0x33536d98, 0xc6ee6618, 0x3c674b25, 0x91ce051e, 0x4a772c94, 0x55e2d70d, 0x70733534, 0x6b47f521, 0xd4a8da98, 0x601ae6f6, 0x0c2adb15, 0x894c3e37, 0x4707e9ae, 0x68c591d9, 0xf2fa9285, 0x840d723c, 0xd829416c, 0x44a00767, 0x26626aec, 0x2a421842, 0xd59de39f, 0xdbb67142, 0x749bdef6, 0x3ec0a6e2, 0x1b192cdf, 0xace59118, 0xdedc3356, 0x1d98fa11, 0xf8b07a47, 0x25edb1e1, 0x80dbda91, 0xe06fe441, 0x8fe7b38e, 0x0eeed70a, 0xa4d621eb, 0xe79c6e5e, 0x7860d8a3, 0x4fa68bfa, 0x545ea673, 0xf58fbf06, 0x25cc1534, 0x87581e6d, 0x5e6adad6, 0xcdc55a31, 0x72e6eb93, 0x2b39b625, 0xb6648db0, 0x8fdeb861, 0x6c216bfa, 0x1e04eb77, 0xec6cedeb, 0x4de13925, 0xb1afeed9, 0x36bd7737, 0x30802277, 0x31cb20e2, 0x59ba341a, 0x84cf4eef, 0xf18c3339, 0x24787414, 0xf693025b, 0xc78134bc, 0x6288009a, 0x58ed793f, 0x55e75dc7, 0x4cea696d, 0x6100a5bc, 0x72c97a23, 0xc21121d6, 0x4484f991, 0xda97a2ba, 0xdbbf6bef, 0xc6853ca7, 0x3d6bdb30, 0xbcf8d2f7,
                                 0x7ca22396, 0xa6bd02f8, 0x7dc6eaa7, 0xbb54c108, 0x03015779, 0x649232dc, 0x02517182, 0xb7814a46, 0x3127912a, 0xd72f4a68, 0x757b7581, 0xeaae5785, 0x6ffbfec7, 0xa5cef31e, 0x48fe345b, 0x8bc45712, 0x8ab258cb, 0xf88a629e, 0x3f9b1d78, 0xc27b4b80, 0x2fdf9d6e, 0x8746683e, 0x143c4586, 0x7ee38888, 0xf9b08521, 0xc2ca48c6, 0xc866427d, 0x44b33b22, 0xaf98d38b, 0x3c4b966c, 0xb3975fb8, 0x00883be5, 0xaae2ab12, 0x060acf39, 0x482ebc51, 0x5ff377f6, 0x7389539a, 0xa5328274, 0x908a7f10, 0x394d97f0, 0x3b0c1b7c, 0x5d4eef8d, 0xf860c81b, 0xa45aff1e, 0x8182a111, 0xba3729ec, 0x97065724, 0xb002776b, 0xf1aeccf1, 0x9a952a1e, 0x63ef11db, 0xb0e9af12, 0xbef8f7ba, 0x6b715dc8, 0xec40ac14, 0xefd75ed4, 0xf3d9e2d2, 0x11d1b57a, 0xc62a6776, 0x0b00e5fb, 0x259f1d4f, 0x6104635d, 0x9bcd2525, 0x79e7dde7, 0x45798c03, 0xb6c20df7, 0xd990342e, 0xa593ddf1, 0x5ef08c09, 0x93aba6e8, 0xe3d6effd, 0xf445745b, 0x11aab09b, 0xc29fd37d, 0xeb67cf9f, 0x4f7dd440, 0xf0783644, 0x9dd60a95, 0xb120ffd1, 0xb1dbaa59,
                                 0x794b35d7, 0x1fc55edf, 0x4ecc75fb, 0x32b79b9b, 0xad7264fe, 0x20456ceb, 0xcc3863e9, 0xb231d06e, 0x6644d208, 0x99651fdb, 0x0a296377, 0x5fdaa60c, 0x59ed3ef8, 0xacd80e57, 0xdd49e6ca, 0x214fffbf, 0x02f138cf, 0xeebca957, 0x7ea37f34, 0xa9ffb051, 0x57674672, 0x87c88b97, 0x4e901b75, 0x13c27062, 0x0bd0a891, 0xbe49de1c, 0xb17b8fec, 0x1981c71c, 0xf81cfa98, 0x76049db0, 0xec334bd9, 0xf157336c, 0xa341c9e5, 0xd96e6a26, 0x58aa4af2, 0xf7d3bd66, 0x27f10b9e, 0x1a338363, 0x6bbc501a, 0x1c730cda, 0x89a5323e, 0xf1e77a82, 0x0be05240, 0xb03956a7, 0x5e6764ee, 0xcdd3a408, 0xefd1fbbf, 0x341aa46d, 0x9bc6e31e, 0x89e4d7b2, 0xcabc6c5b, 0xcac7cb1f, 0x6a80f95b, 0x08aa9d85, 0x9d3fb91d, 0x8b2415e2, 0x3d0964d3, 0xc4ffbc7a, 0x476cc9e3, 0x53f32ffe, 0x59903421, 0x8cafb2b0, 0x8504ea37, 0x5b7f770e, 0xda29a1d5, 0x0b7f0d26, 0x31f8d66d, 0xc1c7463e, 0x92b2708f, 0x1dacc08c, 0x3a9062e7, 0xb5aafe79, 0x58db400d, 0xf2dcf52e, 0x0814272a, 0x9a235a89, 0x7e0f85de, 0x6075c53c, 0xe4aed7d2, 0x716a74ae,
                                 0xaadad3ba, 0xfc830198, 0xa761abe4, 0x3e433916, 0x1c39db46, 0xb2b5cd09, 0x97059521, 0x999e6f10, 0xebb70f89, 0xea6017f1, 0x3e74a51a, 0x7eb47ed8, 0xeb2099e6, 0xd98acf13, 0x0fccc8f6, 0x87e90151, 0xe902acc1, 0x3d816ddf, 0xf35d691b, 0x370bce3c, 0xf9fc978e, 0x3b945fd3, 0xd2bd173d, 0x22c7b209, 0x07aa94fe, 0x7972be71, 0x6afaef2a, 0x1ad0ed3e, 0x2d91b032, 0xaa4ae247, 0x97bbdaeb, 0x3343385e, 0xe45f4679, 0xbcb966e9, 0x28eae152, 0x4ecd77c9, 0x5df35ffe, 0xea0e0f7c, 0xcc849d26, 0x24879aaf, 0xb3d93562, 0x743aa094, 0xeef1ef41, 0x340f76e4, 0xfff0804e, 0x8fcc755f, 0x4fae65ca, 0xbd9d7c2a, 0xacf44590, 0xa63531d8, 0x50beefc4, 0x5665240c, 0x383359e9, 0x785e2fd9, 0x7d2f1b66, 0xb0df64a8, 0xc9c3c474, 0xf78b26d1, 0x9d5036e4, 0x463961b8, 0x633f6226, 0x6388f274, 0x0a2a32ec, 0x4097b6e0, 0x3a409086, 0xc2d0c1fa, 0x982b3c43, 0x374d7303, 0xe89a0533, 0x43ea2350, 0x3915fcf6, 0x3c98735f, 0xa887027b, 0xed2170a5, 0x3b2db962, 0x4effddc8, 0xec6ed59d, 0x4c4690a1, 0x36fd97bc, 0x603de80b,
                                 0xf03853e6, 0xebf3d492, 0xb23dbf25, 0x9cbbb20e, 0x3b8e8647, 0x384d82d2, 0x61bc5f4f, 0x1bd4184b, 0x8780ebe1, 0x3e759c0e, 0xd3c961d9, 0x66e6c2b1, 0x8c394453, 0xb650ce7e, 0x9b2b5061, 0x296af19d, 0x873575b5, 0x6670680c, 0xe27e595c, 0x4591861b, 0x26ec019f, 0xcd82dee6, 0xcc1a994c, 0xcb061fe8, 0x13f66e40, 0x94a9c00a, 0xa8863dba, 0x19e3f04c, 0x9503928e, 0xc4d430e9, 0x4a95ef2d, 0x72865291, 0x1600b44b, 0x3c9d285a, 0x37a530d5, 0x757c099e, 0xa115c984, 0xa8d13fdb, 0x273ef961, 0xc1ed0189, 0xf8ee8457, 0xa61c5e94, 0xc23008af, 0xb96c49a8, 0x739ae779, 0x080680ca, 0xc2b5153a, 0x4b1b9354, 0x124c76ff, 0xc0b724c3, 0x124a9fae, 0xb8a4d236, 0xbda31186, 0xedba3a66, 0xd03c5efa, 0x4104cb96, 0x8ebe3654, 0x81a82cdb, 0x2e202d5b, 0x7194e739, 0x24c9015e, 0x70a941cf, 0x097af8d7, 0x17ddc3f5, 0x5e5994d5, 0x3fff5e60, 0x9e578d6b, 0x59f1b31c, 0x73346480, 0xc64c57a1, 0xd96ff087, 0x1622d287, 0xa6ffa103, 0xdd2806fa, 0x8a911377, 0x278523f9, 0x5e2056c5, 0x511faee6, 0x6ea69dd9, 0x18376132,
                                 0x17594807, 0xaa8186e9, 0xd033c06b, 0xf10b4205, 0x109c43d9, 0x63571f41, 0x7e520442, 0x0d83d15a, 0x9a1bb106, 0x54c7ac4b, 0x429b15b0, 0xac496fcc, 0x5a9ea8d3, 0x18657aa6, 0xbf33dd77, 0x68124bda, 0x13b9f510, 0xf7bfd877, 0xabe9e704, 0x30a33166, 0x815cf11d, 0x39175c95, 0x18756d44, 0x1dc8256d, 0x726ef86b, 0x2f2e181e, 0x96f864ca, 0x9c5a544c, 0x5e7cf8fe, 0x71078074, 0xd6502936, 0x0809cfa0, 0x6806c6d5, 0xc970d239, 0x16e27ea5, 0x5180589b, 0x00cc961b, 0x018a93fa, 0x81e864bf, 0x8cc4cfce, 0xee6704e1, 0x29ac80b1, 0x9d6f1f23, 0x944bda80, 0x22391ab0, 0x4cad23ec, 0x7300a7bf, 0xe6a52bc0, 0x696df14d, 0x6b564674, 0x24f9e72f, 0x6d8717a9, 0x50a9a041, 0xabfbd168, 0x527b36b7, 0xc60b338c, 0x4660a6d2, 0x84648db2, 0xedf03551, 0xd27185b2, 0xc1c3f1e1, 0x9b7ba492, 0x643f52c6, 0xe6f04539, 0xccc4b5c9, 0x1fb22cd6, 0x5a976ece, 0x54b69dcb, 0xd5480451, 0x1b994bac, 0xa347066e, 0xa66666be, 0xe3032f68, 0xf53d88c9, 0xfaed2b34, 0x33185648, 0x0bccd9a7, 0x076ba547, 0x61aa62f6, 0x7c070530,
                                 0x68d859ed, 0x68c11e92, 0xcf40acc1, 0xcffe2f23, 0xf74e366e, 0x06916db0, 0xd187d41d, 0xa602d58f, 0x3dc2f9bc, 0xdc9f9ed2, 0x53e58c31, 0xa39fe4d3, 0x4b7e84d4, 0x361f8951, 0xbcd628a6, 0xde7616a9, 0xab651e05, 0x65501d64, 0x3f30f50f, 0xf5987638, 0xec2f92e2, 0xa64323b6, 0x3fc2c752, 0x69991563, 0x9b7eba8f, 0x840691a9, 0xa6edb785, 0x10f33201, 0xb711f99d, 0x2e10a9d4, 0xe41531ca, 0xf8dc60a0, 0x07572d8e, 0xb003d117, 0x7f1214ca, 0x7b4bada4, 0x52bf12c6, 0x8a0b2eff, 0x5470e9bc, 0x2bcfdfa2, 0xa2378880, 0xb74ae264, 0x99a6f4c9, 0xc3ac3530, 0xdc6d9764, 0x782b1123, 0x3ff39bc0, 0x49d6f4cf, 0x58dc76db, 0x6ccda92e, 0x1a2e0a6e, 0x5e72d59a, 0x4228715a, 0x52f376d3, 0xcdb82679, 0x3fc75512, 0xcf11e397, 0xb022a12d, 0xe96cedab, 0xd2ed52aa, 0x7e7c9c6b, 0xb4836928, 0xfab50715, 0xdd5f5efa, 0xb44a597a, 0x6217e532, 0x5b42f348, 0x6d3221fb, 0x6a6b1ea1, 0x56dd495b, 0xb77dc9d1, 0xf87c6a2b, 0xd1c37d65, 0x6fbf8130, 0x3b605e2d, 0xef5472f8, 0xdb4e64b1, 0x2ae5f61e, 0xee5c8a7e, 0x15e52493,
                                 0x489740c5, 0x752da351, 0xebddfb31, 0x4c3b7583, 0xf5179eed, 0xfc5b04ed, 0x75a968ff, 0xea05c892, 0xa2bd8879, 0x9e250008, 0x91d9f98d, 0xd381a7b8, 0x5bdd3927, 0xda1aea0c, 0xbd05b57a, 0x7cb92f7f, 0x23cacb1b, 0xe2c36bd3, 0x88d8e1a8, 0x5a3c3db6, 0x5748987c, 0x3d8b0b0c, 0x8f7383ff, 0xff79d80c, 0x9ee75f81, 0x3b6ab3e5, 0xd9331a19, 0x4c074f69, 0xec3de15c, 0x219a0a4d, 0xcacab500, 0x2991176e, 0x724b33f8, 0x7bb3b110, 0xd9f2a6b4, 0x78c8e4cb, 0x09f41829, 0xd3603000, 0x4aeb73ad, 0x69bf4dca, 0x53d6d606, 0xaa43dbc1, 0x78c3811e, 0x19149496, 0xcd78ec77, 0x4161e9a7, 0x2fa48ba9, 0xdc69dbe7, 0x238ed0ad, 0xc8b39fc0, 0x69fe8967, 0x1b32d92f, 0xc27bf5c3, 0xc89ef38a, 0x44fc7662, 0x6a975472, 0xb5089631, 0x22b448c5, 0x3ffed38e, 0xf53b201a, 0xed517253, 0x54526313, 0x2d6f8654, 0x55457c3c, 0xbd81b4e5, 0x703e30e8, 0x81c28366, 0x7ecd682f, 0x7023c7eb, 0x762f6f2f, 0x5baf9e59, 0x481aba14, 0xa52d42d5, 0x4c2d353b, 0x940fae85, 0x67491e55, 0xd4c4b2c9, 0xecabf54d, 0x4853ae97};

static const unsigned bufferSize1279 = 1279;
static const unsigned tauswortheQ1279 = 1063;

RandomKernel1279::RandomKernel1279(unsigned seed, Verbosity verb) : Tausworthe(bufferSize1279, tauswortheQ1279, seedbuf1279, seed, verb) {}

RandomKernel1279::RandomKernel1279(const char *locIdent, unsigned counterOverride, Verbosity verb) : Tausworthe(bufferSize1279, tauswortheQ1279, seedbuf1279, locIdent, counterOverride, verb) {}