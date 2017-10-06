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
    public class NetworkReconstructionProcessing : PluginInterop.Python.NetworkFromMatrix
    {
        public override string Heading => "Reconstruct";
        public override string Name => "ANAT";
        public override string Description => "Reconstruct a signaling pathway";
        public override Bitmap2 DisplayImage => GraphUtils.ToBitmap2(Properties.Resources.icon);

        protected override string[] ReqiredPythonPackages => new[] { "perseuspy", "phos" };


        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            byte[] code = (byte[])Properties.Resources.ResourceManager.GetObject("reconstruct_network");
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
                new SingleChoiceWithSubParams("Select terminals from")
                {
                    Value = Math.Max(0, mdata.CategoryColumnNames.FindIndex(col => col.ToLower().Contains("significant"))),
                    Values = mdata.CategoryColumnNames,
                    Help = "Select the end points (terminals) of the signaling networks.",
                    SubParams = Enumerable.Range(0, mdata.CategoryColumnCount).Select(i =>
                    {
                        var subParam = new MultiChoiceParam("Select", new int[0])
                        {
                            Values = mdata.GetCategoryColumnValuesAt(i)
                        };
                        return new Parameters(subParam);
                    }).ToArray()
                }, 
                new FileParam("Network", networkLocation) {Help="Weighted PPI network."},
                new DoubleParam("Network confidence", 0.5) {Help = "Imposes a confidence cutoff on the edges of the PPI network."},
                new IntParam("Network node degree threshold", 150) {Help = "Removes highly connected nodes (such as Ubiquitin) from the network."},
            };
        }
    }
}