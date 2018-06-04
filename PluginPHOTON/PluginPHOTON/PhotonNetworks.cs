using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using BaseLib.Graphic;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Matrix;
using Path = System.IO.Path;

namespace PluginPHOTON
{
    public class PhotonNetworks : PluginInterop.Python.NetworkFromMatrix
    {
        public override string Heading => "Modifications";
        public override string Name => "PHOTON";
        public override string Description => "Reconstruct a signaling pathway";
        public override Bitmap2 DisplayImage => GraphUtils.ToBitmap2(Properties.Resources.icon);

        public override int NumSupplTables => 1;
        public override string HelpOutput => "Reconstructed signaling networks";
        public override string[] HelpSupplTables => new [] {"Signaling score for all proteins in the PPI network with sufficient data."};
        protected override string[] ReqiredPythonPackages => new[] { "perseuspy", "phos" };


        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            byte[] code = (byte[])Properties.Resources.ResourceManager.GetObject("pipeline");
            codeFile = Path.GetTempFileName();
            File.WriteAllText(codeFile, Encoding.UTF8.GetString(code));
            return true;
        }

        protected override Parameter[] SpecificParameters(IMatrixData mdata, ref string errString)
        {
            var networkLocation = string.Empty;
            try
            {
                networkLocation = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "conf", "H_sapiens.net");
            }
            catch (Exception ex)
            {
                Debug.WriteLine(ex.Message);
            }
            return new Parameter[]
            {
                new SingleChoiceParam("GeneID")
                {
                    Value = Math.Max(0, mdata.StringColumnNames.FindIndex(col => col.ToLower().Equals("geneid"))),
                    Values = mdata.StringColumnNames,
                    Help = "Single entrez gene id per row"
                },
                new StringParam("Signaling source")
                {
                    Help = "Has to be a human entrez gene id. Select the starting point of the signaling network (optional)."
                }, 
                new SingleChoiceParam("Amino acid")
                {
                    Value = mdata.CategoryColumnNames.FindIndex(col => col.ToLower().Equals("amino acid")),
                    Values = mdata.CategoryColumnNames,
                    Help = "Amino acid of the modification"
                }, 
                new SingleChoiceParam("Position")
                {
                    Value = mdata.NumericColumnNames.FindIndex(col => col.ToLower().Equals("position")),
                    Values = mdata.NumericColumnNames,
                    Help = "Position of the modification within the protein"
                },
                new FileParam("Network", networkLocation) {Help="Weighted PPI network."},
                new DoubleParam("Network confidence", 0.5) {Help = "Imposes a confidence cutoff on the edges of the PPI network."},
                new IntParam("Network node degree threshold", 150) {Help = "Removes highly connected nodes (such as Ubiquitin) from the network."},
                new IntParam("Required number of observations", 4) {Help = "Required minumum number of observations for score calculation."},
                new IntParam("Number of permutations", 1000) {Help = "Number of permutations used for empirical p-value calculation"},
            };
        }
    }
}