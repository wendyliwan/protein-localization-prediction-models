import matplotlib.pyplot as plt
from PIL import Image
import torch
from PIL import ImageFont, ImageDraw

if __name__ == "__main__":
    # 有 GPU 就用 GPU，没有就用 CPU
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print('device', device)

    # 导入中文字体，指定字体大小
    #font = ImageFont.truetype('SimHei.ttf', 50)

    from torchvision.models import resnet18
    model = resnet18(pretrained=True).eval().to(device)

    from torchcam.methods import SmoothGradCAMpp
    # CAM GradCAM GradCAMpp ISCAM LayerCAM SSCAM ScoreCAM SmoothGradCAMpp XGradCAM
    cam_extractor = SmoothGradCAMpp(model)

    from torchvision import transforms
    # 测试集图像预处理-RCTN：缩放、裁剪、转 Tensor、归一化
    test_transform = transforms.Compose([transforms.Resize(512),
                                     transforms.CenterCrop(224),
                                     transforms.ToTensor(),
                                     transforms.Normalize(
                                         mean=[0.485, 0.456, 0.406],
                                         std=[0.229, 0.224, 0.225])
                                    ])
    img_path = 'E:\Testis_IHC\Breast/cancer_all\CCNT1/CCNT1_2392_Duct carcinoma_.png'
    img_pil = Image.open(img_path)
    input_tensor = test_transform(img_pil).unsqueeze(0).to(device) # 预处理
    pred_logits = model(input_tensor)
    pred_top1 = torch.topk(pred_logits, 1)
    pred_id = pred_top1[1].detach().cpu().numpy().squeeze().item()
    activation_map = cam_extractor(pred_id, pred_logits)
    activation_map = activation_map[0][0].detach().cpu().numpy()
    plt.imshow(activation_map)
    #plt.show()
    from torchcam.utils import overlay_mask
    result = overlay_mask(img_pil, Image.fromarray(activation_map), alpha=0.7)
    plt.imshow(result)
    plt.show()